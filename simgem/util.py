"""Logging, subprocesses and locating the compiled Garfield programs."""
import contextlib
import datetime
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class _Tee:
    """Write to several streams at once (console + a log file)."""

    def __init__(self, *streams):
        self._s = streams

    def write(self, x):
        for s in self._s:
            s.write(x)
            s.flush()

    def flush(self):
        for s in self._s:
            s.flush()

    def fileno(self):
        return self._s[0].fileno()

    def isatty(self):
        return False


@contextlib.contextmanager
def stage_log(rundir, name):
    """Tee stdout/stderr of a stage into <rundir>/<name>.log."""
    os.makedirs(rundir, exist_ok=True)
    path = os.path.join(rundir, f"{name}.log")
    f = open(path, "w")
    old_out, old_err = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = _Tee(old_out, f), _Tee(old_err, f)
    print(f"===== {name}  {datetime.datetime.now().isoformat(timespec='seconds')} =====")
    try:
        yield path
    finally:
        sys.stdout, sys.stderr = old_out, old_err
        f.close()
        old_out.write(f"[log -> {path}]\n")


def run_streamed(cmd, cwd=None, log_path=None, echo=True, keep=None):
    """Run a subprocess. Its output goes to `log_path` (if given) and, when
    `echo`, to stdout; `keep(line)` can filter what is shown/logged.
    Raises CalledProcessError (with the log path) on failure."""
    cmd = [str(c) for c in cmd]
    print("+ " + " ".join(cmd) + (f"   [log: {log_path}]" if log_path else ""))
    lf = open(log_path, "w") if log_path else None
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True, bufsize=1)
    tail = []
    for line in p.stdout:
        if keep and not keep(line):
            continue
        if lf:
            lf.write(line)
        if echo:
            sys.stdout.write(line)
        tail = (tail + [line])[-30:]
    p.wait()
    if lf:
        lf.close()
    if p.returncode != 0:
        if not echo:
            sys.stdout.write("".join(tail))
        raise subprocess.CalledProcessError(p.returncode, cmd)


def binary(name):
    """Path of a compiled SimGEM program: $SIMGEM_BIN, ./build, then PATH."""
    cands = []
    if os.environ.get("SIMGEM_BIN"):
        cands.append(os.path.join(os.environ["SIMGEM_BIN"], name))
    cands.append(os.path.join(ROOT, "build", name))
    for c in cands:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    w = shutil.which(name)
    if w:
        return w
    sys.exit(f"Cannot find '{name}'. Build it first:\n"
             f"  cmake -S . -B build && cmake --build build -j")


def _condor_cpus():
    """CPUs of the HTCondor slot (grid job), or None outside HTCondor. HTCondor
    sets OMP_NUM_THREADS to request_cpus; the job ad is the fallback."""
    if not (os.environ.get("_CONDOR_JOB_AD") or os.environ.get("_CONDOR_SCRATCH_DIR")):
        return None
    if os.environ.get("OMP_NUM_THREADS", "").isdigit():
        return int(os.environ["OMP_NUM_THREADS"])
    try:
        with open(os.environ["_CONDOR_JOB_AD"]) as f:
            for line in f:
                k, _, v = line.partition("=")
                if k.strip().lower() == "requestcpus" and v.strip().isdigit():
                    return int(v.strip())
    except (KeyError, OSError):
        pass
    return 1


def ncpu():
    """Default number of parallel workers: $SIMGEM_JOBS, else the CPUs of the
    HTCondor slot, else the cores this process may use."""
    env = os.environ.get("SIMGEM_JOBS")
    if env:
        return max(1, int(env))
    n = _condor_cpus()
    if n:
        return max(1, n)
    if hasattr(os, "sched_getaffinity"):
        return len(os.sched_getaffinity(0)) or 1
    return os.cpu_count() or 4
