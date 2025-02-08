import argparse, os, shutil

def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input file with line with ids", required=True)
    parser.add_argument('--folder', help="folder with files", required=True)
    parser.add_argument("-o", "--output", help="output folder with ", required=True)
    args = parser.parse_args()
    return args

def find_file(folder, query)->list:
    for file in os.listdir(folder):
        if file.split('.')[0] == query:
            return f"{folder}/{file}"
    return None

def main():
    args = arguments()
    with open(args.input, 'r') as f:
        for line in f:
            line = line.strip()
            file_path = find_file(args.folder, line)
            output_file_path = f"{args.output}/{file_path.split('/')[-1]}"

            os.makedirs(args.output, exist_ok=True)
            shutil.copy(file_path, output_file_path)

main()