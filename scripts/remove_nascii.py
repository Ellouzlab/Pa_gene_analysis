import os
import argparse

def replace_non_ascii(text):
    return ''.join(char if ord(char) < 128 else '?' for char in text)

def process_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        
        if os.path.isfile(input_path):
            with open(input_path, 'r', encoding='utf-8', errors='ignore') as infile:
                content = infile.read()
            corrected_content = replace_non_ascii(content)
            with open(output_path, 'w', encoding='utf-8') as outfile:
                outfile.write(corrected_content)
            print(f'Processed {filename}')

def main():
    parser = argparse.ArgumentParser(description='Replace non-ASCII characters in files')
    parser.add_argument('--input', required=True, help='Input directory containing files')
    parser.add_argument('--output', required=True, help='Output directory for corrected files')
    
    args = parser.parse_args()
    
    process_files(args.input, args.output)

if __name__ == '__main__':
    main()
