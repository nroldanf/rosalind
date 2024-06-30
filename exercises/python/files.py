if __name__ == "__main__":
    # Given: A file containing at most 1000 lines.
    # Return: A file containing all the even-numbered lines from the original file. Assume 1-based numbering of lines.
    with open("../../data/rosalind_ini5.txt", "r") as file:
        lines = file.readlines()
    with open("../../data/rosalind_ini5_output.txt", "w") as output_file:
        for i, line in enumerate(lines, start=1):
            if i % 2 == 0:
                output_file.write(line)