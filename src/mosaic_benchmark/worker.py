"""Single-threaded resource wrapper: avoid preexec_fn inside scheduler threads."""
import argparse
import os
import resource


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cpus', required=True)
    p.add_argument('--memory-bytes', required=True, type=int)
    p.add_argument('--output-bytes', required=True, type=int)
    p.add_argument('--stack-bytes', type=int)
    p.add_argument('command', nargs=argparse.REMAINDER)
    args = p.parse_args()
    command = args.command[1:] if args.command and args.command[0] == '--' else args.command
    if not command:
        p.error('Native command missing')
    os.sched_setaffinity(0, [int(x) for x in args.cpus.split(',')])
    resource.setrlimit(resource.RLIMIT_AS, (args.memory_bytes, args.memory_bytes))
    resource.setrlimit(resource.RLIMIT_FSIZE, (args.output_bytes, args.output_bytes))
    if args.stack_bytes:
        resource.setrlimit(resource.RLIMIT_STACK, (args.stack_bytes, args.stack_bytes))
    os.execv(command[0], command)


if __name__ == '__main__':
    main()
