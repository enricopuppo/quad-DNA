#!/usr/bin/env python3

import argparse
import datetime as dt
import errno
import io
import os
import pty
import re
import signal
import subprocess
import sys
import traceback
from pathlib import Path


def find_input_objs(data_dir: Path) -> list[Path]:
	"""Find OBJ files where directory name matches file name (e.g. XXX/XXX.obj)."""
	matches: list[Path] = []
	for dirpath, _, _ in __import__("os").walk(data_dir):
		current_dir = Path(dirpath)
		candidate = current_dir / f"{current_dir.name}.obj"
		if candidate.is_file():
			matches.append(candidate)
	matches.sort()
	return matches


def resolve_executable(executable_arg: str) -> Path:
	exe_path = Path(executable_arg)
	if exe_path.exists():
		return exe_path.resolve()

	# Fall back to PATH lookup if the provided value is a command name.
	import shutil

	found = shutil.which(executable_arg)
	if found is None:
		raise FileNotFoundError(
			f"Executable not found: '{executable_arg}'. "
			"Pass --exe with an absolute/relative path or ensure it is in PATH."
		)
	return Path(found).resolve()


def write_separator(logf, title: str) -> None:
	bar = "=" * 100
	logf.write(f"\n{bar}\n")
	logf.write(f"{title}\n")
	logf.write(f"{bar}\n")


def run_with_pty(command: list[str]) -> tuple[int, str]:
	"""Run a command attached to a pseudo-terminal to capture console-emitted diagnostics."""
	master_fd, slave_fd = pty.openpty()
	try:
		proc = subprocess.Popen(
			command,
			stdin=subprocess.DEVNULL,
			stdout=slave_fd,
			stderr=slave_fd,
			close_fds=True,
		)
	finally:
		os.close(slave_fd)

	chunks: list[bytes] = []
	while True:
		try:
			data = os.read(master_fd, 4096)
			if not data:
				break
			chunks.append(data)
		except OSError as exc:
			# On macOS, EIO from PTY read means the child has closed the terminal.
			if exc.errno == errno.EIO:
				break
			raise

	os.close(master_fd)
	returncode = proc.wait()
	output = b"".join(chunks).decode("utf-8", errors="replace")
	return returncode, output


def extract_failure_lines_from_side_log(input_file: Path) -> str | None:
	"""Return only failure-relevant lines from <input>.log, if present."""
	log_path = input_file.with_suffix(input_file.suffix + ".log")
	if not log_path.is_file():
		return None

	try:
		content = log_path.read_text(encoding="utf-8", errors="replace")
	except Exception:
		return None

	all_lines = [line for line in content.splitlines() if line.strip()]
	if not all_lines:
		return None

	pattern = re.compile(r"assert|abort|error|failed|failure|fatal", re.IGNORECASE)
	filtered = [line for line in all_lines if pattern.search(line)]
	if not filtered:
		filtered = all_lines[-3:]

	return "\n".join(filtered)


def run_batch(data_dir: Path, exe_path: Path, log_path: Path, only_failures_log: bool = False) -> int:
	input_files = find_input_objs(data_dir)

	with log_path.open("w", encoding="utf-8") as logf:
		logf.write(f"DNA batch run: {dt.datetime.now().isoformat()}\n")
		logf.write(f"Data directory: {data_dir}\n")
		logf.write(f"Executable: {exe_path}\n")
		logf.write(f"Found {len(input_files)} input files\n")
		logf.write(f"Only failures in per-file log: {only_failures_log}\n")

		if not input_files:
			write_separator(logf, "NO INPUT FILES FOUND")
			logf.write(
				"No files matched the expected pattern <subdir>/<subdir>.obj under the data directory.\n"
			)
			return 0

		total_inputs = len(input_files)
		failures = 0
		logged_sections = 0
		for idx, input_file in enumerate(input_files, start=1):
			header = f"[{idx}/{total_inputs}] INPUT: {input_file}"
			entry = io.StringIO()
			write_separator(entry, header)

			try:
				returncode, console_output = run_with_pty([str(exe_path), str(input_file)])
			except Exception as exc:
				failures += 1
				entry.write(f"ERROR launching process: {exc}\n")
				entry.write("Traceback:\n")
				entry.write(traceback.format_exc())
				if not traceback.format_exc().endswith("\n"):
					entry.write("\n")
				logf.write(entry.getvalue())
				logged_sections += 1
				print(f"[{idx}/{total_inputs}]", flush=True)
				continue

			entry.write(f"Exit code: {returncode}\n")
			entry.write("Program stdout (console-captured):\n")
			if console_output:
				entry.write(console_output)
				if not console_output.endswith("\n"):
					entry.write("\n")
			else:
				entry.write("<no output>\n")

			if returncode != 0:
				if not console_output or not console_output.strip():
					entry.write("Failure diagnostics:\n")
					side_log_diag = extract_failure_lines_from_side_log(input_file)
					if side_log_diag:
						entry.write("Relevant lines from program side log:\n")
						entry.write(side_log_diag)
						if not side_log_diag.endswith("\n"):
							entry.write("\n")
					if returncode < 0:
						sig = -returncode
						try:
							sig_name = signal.Signals(sig).name
						except ValueError:
							sig_name = "UNKNOWN"
						entry.write(f"Process terminated by signal {sig} ({sig_name}).\n")
					elif returncode >= 128:
						sig = returncode - 128
						try:
							sig_name = signal.Signals(sig).name
						except ValueError:
							sig_name = "UNKNOWN"
						entry.write(
							f"Shell reported abnormal termination, likely signal {sig} ({sig_name}).\n"
						)
					else:
						entry.write(
							"Process exited with non-zero code but produced no console output.\n"
						)
				failures += 1

			if (not only_failures_log) or (returncode != 0):
				logf.write(entry.getvalue())
				logged_sections += 1

			print(f"[{idx}/{total_inputs}]", flush=True)

		write_separator(logf, "BATCH SUMMARY")
		logf.write(f"Total inputs: {len(input_files)}\n")
		logf.write(f"Logged input sections: {logged_sections}\n")
		logf.write(f"Failures: {failures}\n")
		logf.write(f"Succeeded: {len(input_files) - failures}\n")

	return 1 if failures else 0


def main() -> int:
	parser = argparse.ArgumentParser(
		description=(
			"Run dna_quad on all matching OBJ inputs in a dataset directory tree. "
			"Expected pattern: <subdir>/<subdir>.obj"
		)
	)
	parser.add_argument("data_dir", help="Path to the dataset root directory")
	parser.add_argument(
		"--exe",
		default="build/dna_quad",
		help="Path or command name for the dna_quad executable (default: build/dna_quad)",
	)
	parser.add_argument(
		"--log",
		default=None,
		help="Path to output log file (default: <data_dir>/DNA_batch.log)",
	)
	parser.add_argument(
		"--only-failures-log",
		action="store_true",
		help="Write per-file sections only for inputs with non-zero exit code",
	)
	args = parser.parse_args()

	data_dir = Path(args.data_dir).expanduser().resolve()
	if not data_dir.is_dir():
		print(f"Error: data_dir is not a directory: {data_dir}", file=sys.stderr)
		return 2

	try:
		exe_path = resolve_executable(args.exe)
	except FileNotFoundError as exc:
		print(f"Error: {exc}", file=sys.stderr)
		return 2

	log_path = Path(args.log).expanduser().resolve() if args.log else data_dir / "DNA_batch.log"
	log_path.parent.mkdir(parents=True, exist_ok=True)

	exit_code = run_batch(
		data_dir=data_dir,
		exe_path=exe_path,
		log_path=log_path,
		only_failures_log=args.only_failures_log,
	)
	print(f"Batch completed. Log written to: {log_path}")
	return exit_code


if __name__ == "__main__":
	raise SystemExit(main())
