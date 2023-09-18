from os.path import join, splitext
import sys


def split_file(file_path, n):
    file = open(file_path, 'r')
    lines = file.readlines()
    num_lines = len(lines)
    num_lines_per_file = int(num_lines / n)

    for i in range(n):
        fileOutputName = splitext(file_path)[0] + '_part' + str(i) + '.csv'
        fileOutput = open(fileOutputName, 'w')
        start = i * num_lines_per_file
        if i == n - 1:
            end = num_lines + 1
        else:
            end = start + num_lines_per_file
        for line in lines[start:end]:
            fileOutput.write(line)
        fileOutput.close()

    file.close()


if __name__ == '__main__':
    split_file(sys.argv[1], int(sys.argv[2]))
