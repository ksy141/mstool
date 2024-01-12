import argparse
import mstool
import sqlite3

description = "dump a DMS file into a SQL file"

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Positional arguments
    parser.add_argument("input",  help="provide a DMS file")
    parser.add_argument("output", help="provide a name of a SQL file")

    return parser.parse_args()

def main():
    args = parse_args()
    fmt  = args.input.split('.')[-1]
    assert fmt.lower() == 'dms', 'Should provide a DMS file'

    conn = sqlite3.connect(args.input)
    with open(args.output, 'w') as f:
        for line in conn.iterdump():
            f.write('%s\n' %line)

if __name__ == '__main__':
    main()

