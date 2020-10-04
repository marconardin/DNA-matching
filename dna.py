from sys import argv, exit
import csv


def main(argv):
    if len(argv) != 3:
        print("Usage: python dna.py data.csv sequence.txt")
        exit(1)

    csv_file = argv[1]
    sequence_file = argv[2]

    with open(csv_file, 'r') as csvf:
        reader = csv.reader(csvf)
        dna_strs = next(reader)[1:]
        entries = [line for line in reader]

    # create dictionary of subject and DNA STR counts
    csv_dict = matrix_to_dict(dna_strs, entries)

    with open(sequence_file, 'r') as seqf:
        sequence = seqf.read().rstrip()

    # find repeats
    str_counts = get_str_repeats(dna_strs, sequence)

    # find potential match to subject
    match = match_repeats_to_person(str_counts, csv_dict)

    print(match)


def matrix_to_dict(dna_strs, lines):
    names = [line[0] for line in lines]

    _dict = {}
    for i, name in enumerate(names):
        for j, str_seq in enumerate(dna_strs):
            if name not in _dict:
                _dict[name] = {}
            _dict[name][str_seq] = int(lines[i][j+1])

    return _dict


def get_str_repeats(dna_strs, main_sequence):
    str_counts = {}

    for strseq in dna_strs:
        strseq_len = len(strseq)
        strseq_max = 0
        strseq_curr_count = 0

        i = 0
        while i < len(main_sequence):
            # check if current substring (length of STR) matches the DNA sequence at this index
            if main_sequence[i:i+strseq_len] == strseq:
                strseq_curr_count += 1
                if strseq_curr_count > strseq_max:
                    strseq_max = strseq_curr_count
                # skip ahead the length of the STR
                i += strseq_len
            else: # STR not found, reset the count and go to next base
                strseq_curr_count = 0
                i += 1

        str_counts[strseq] = strseq_max

    return str_counts


def match_repeats_to_person(str_counts, csv_dict):
    for subject, subject_counts in csv_dict.items():
        if str_counts == subject_counts:
            return subject
    return "No match"

if __name__ == "__main__":
    main(argv)