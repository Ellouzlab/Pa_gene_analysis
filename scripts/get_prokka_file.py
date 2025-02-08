import argparse, shutil, os

def args():
    '''
    Get arguments
    :return:
    '''
    parser = argparse.ArgumentParser(description = 'Get prokka files from a directory if they have a specific extension')
    parser.add_argument('-i', '--input', required = True, help = 'Input directory')
    parser.add_argument('-o', '--output', required = True, help = 'Output directory')
    parser.add_argument('-e', '--extension', required = True, help = 'File extension')
    return parser.parse_args()

def get_prokka_file():
    '''
    Get prokka files from a directory and copy them to a new directory if they have a specific extension
    :return:
    '''
    arguments = args()

    if not os.path.exists(arguments.output):
        os.mkdir(arguments.output)

    for folder in os.listdir(arguments.input):
        folder_path = os.path.join(arguments.input, folder)
        for file in os.listdir(folder_path):
            if file.endswith(arguments.extension):
                file_path = os.path.join(folder_path, file)
                new_file_path = f"{arguments.output}/{folder.split('.')[0]}.{arguments.extension}"
                shutil.copy(file_path, new_file_path)

if __name__ == '__main__':
    get_prokka_file()
