import json
import sys

if __name__ == "__main__":

    with open(sys.argv[1]) as data_file:
        conf = json.load(data_file)

    print conf[sys.argv[2]]
