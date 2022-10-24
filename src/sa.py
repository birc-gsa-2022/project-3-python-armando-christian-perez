import argparse
import sys

argparser = argparse.ArgumentParser(
    description="Exact matching using suffix tree construction")
argparser.add_argument("genome", type=argparse.FileType('r'))
argparser.add_argument("reads", type=argparse.FileType('r'))
args = argparser.parse_args()


def fasta_translator(file):
    output_dict = {}
    start = True
    for i in file:
        print("hello", "\t", "i")
        i = i.strip()
        if len(i) == 0:
            continue
        if i[0] == ">":
            if  not start:
                output_dict[name] = seq
            name = i[1:].strip()
            seq = ""
            if start:
                start = False
        elif start:
            continue
        else:
            seq += i
    output_dict[name] = seq
    return output_dict

def fastq_translator(file):
    output_dict = {}
    for i in file:
        if i[0] == "@":
            name = i[1:].strip()
        else:
            seq = i.strip()
            output_dict[name] = seq
    return(output_dict)

def annotate_string(string): # Should be run with a sentinel attached
    LS_array = ["$" for i in range(len(string))]
    LS_array[-1] = "S"

    for i in range(len(string) - 2, -1, -1):
        if string[i] > string[i + 1]:
            LS_array[i] = "L"
        else:
            LS_array[i] = "S"
    return LS_array

def get_LMS(LS_array):
    LMS_array = []

    for i in range(1, len(LS_array)):
        if LS_array[i] == "S" and LS_array[i - 1] == "L":
            LMS_array.append(i)
    return LMS_array

def get_LMS_blocks(string, LMS_array):
    LMS_blocks = []
    LMS_array.append(len(string))
    for i in range(len(LMS_array) - 1):
        LMS_blocks.append(string[LMS_array[i] : LMS_array[i + 1] + 1])
    LMS_array.pop()
    return LMS_blocks

def impute_buckets(string):
    buckets = {}
    for i in range(len(string)):
        if string[i] in buckets:
            buckets[string[i]].append(i)
        else:
            buckets[string[i]] = [i]
    buckets = dict(sorted(buckets.items()))
    return buckets

def output_bucket(buckets):
    bucket_out = {}
    for i in buckets:
        bucket_out[i] = [-1 for i in range(len(buckets[i]))]
    return bucket_out

def first_last_index(bucket, fl):
    if fl == "last":
        for i in range(len(bucket) - 1, -1, -1):
            if bucket[i] == -1:
                return i
    elif fl == "first":
        for i in range(len(bucket)):
            if bucket[i] == -1:
                return i

def first_pass(buckets, bucket_out, LMS_array, LMS_blocks):
    for i in range(len(LMS_array) -1, -1, -1):
        index = first_last_index(bucket_out[LMS_blocks[i][0]], "last")
        bucket_out[LMS_blocks[i][0]][index] = LMS_array[i]
    return bucket_out

def second_pass(bucket_out, LS_array, string):
    for key in bucket_out:
        for element in bucket_out[key]:
            if element < 1:
                continue
            if LS_array[element - 1] == "L":
                index = first_last_index(bucket_out[string[element - 1]], "first")
                bucket_out[string[element - 1]][index] = element - 1
    return bucket_out

def bucket_indices(bucket_out):
    indices = {}
    for key in bucket_out:
        indices[key] = len(bucket_out[key]) - 1
    return indices

def third_pass(bucket_out, LS_array, string):
    indices = bucket_indices(bucket_out)

    for key in sorted(bucket_out, reverse=True):
        for i in range(len(bucket_out[key]) - 1, -1, -1):
            if bucket_out[key][i] < 1:
                continue
            elif LS_array[bucket_out[key][i] - 1] == "S":
                buck = string[bucket_out[key][i] - 1]
                bucket_out[buck][indices[buck]] = bucket_out[key][i] - 1
                indices[buck] -= 1
    return bucket_out

def number_LMS(bucket_out, LMS_blocks, LMS_array):
    reduced = [0 for i in range(len(LMS_array))]
    previous = None
    count = -1
    for key in bucket_out:
        for element in bucket_out[key]:
            if element in LMS_array:
                block = LMS_blocks[LMS_array.index(element)]
                if block != previous:
                    count += 1
                previous = block
                reduced[LMS_array.index(element)] = str(count)
    return "".join(reduced)

def compute_t1(string):
    LS_array = annotate_string(string)
    LMS_array = get_LMS(LS_array)
    LMS_blocks = get_LMS_blocks(string, LMS_array)
    buckets = impute_buckets(string)
    bucket_out = output_bucket(buckets)
    bucket_out = first_pass(buckets, bucket_out, LMS_array, LMS_blocks)
    bucket_out = second_pass(bucket_out, LS_array, string)
    bucket_out = third_pass(bucket_out, LS_array, string)
    reduced = 0#number_LMS(bucket_out, LMS_blocks, LMS_array)
    return reduced, bucket_out

#def SAIS(string):
#    stack = list()
#    while True:
#        reduced, bucket_out = compute_t1
#        if reduced[-1] != len()

#def mostly_naive(string):
#    forget, bucket_out = compute_t1(string)
#    SA = []
#    for i in bucket_out:
#        SA.extend(bucket_out[i])
#    SA = [string[SA[int(i)]:] for i in SA]
#    SA.sort()
#    return(SA)

def completely_naive(string):
    SA = [i for i in range(len(string))]
    SA = [string[SA[int(i)]:] for i in SA]
    SA.sort()
    return(SA)

def test_SA(SA):
    for i in SA:
        print(i)

def LCP_naive(SA, string):
    LCP = [0, 0]
    for i in range(2, len(SA)):
        count = 0
        for j in range(len(string)):
            if string[SA[i] + j] == string[SA[i - 1] + j]:
                count += 1
            else:
                break
        LCP.append(count)
    return LCP

def naive(string):
    string += "$"
    SA = completely_naive(string)
    SA = [len(string) - len(i) for i in SA]
    LCP = LCP_naive(SA, string)
    return SA, LCP

def compare_strings(pattern, suffix):
    loops = min(len(pattern), len(suffix))

    for i in range(loops):
        if pattern[i] < suffix[i]:
            return False
        if pattern[i] > suffix[i]:
            return True

def SA_search(SA, string, pattern):
    lo = 0
    hi = len(SA) - 1
    while True:
        mid = int((lo + hi)/2)
        suffix = string[SA[mid]:]
        match = compare_strings(pattern, suffix)
        if match is None:
            return mid
        elif hi == lo:
            return None
        elif match:
            lo = mid + 1
        elif not match:
            hi = mid - 1

def SA_match(SA, LCP, string, pattern):
    string += "$"
    match = SA_search(SA, string, pattern)
    l = len(pattern)
    if match is None:
        return []
    lo = match
    hi = match + 1

    while hi < len(LCP):
        if LCP[hi] >= l:
            hi += 1
        else:
            break
    while lo > 0:
        if LCP[lo] >= l:
            lo -= 1
        else:
            break
    return SA[lo:hi]

def matches_to_SAM(read_file, reference_file):
    read_name = []
    reference_name = []
    match_index = []
    CIGARS = []
    match_string = []

    fasta_dict, fastq_dict = fasta_translator(reference_file), fastq_translator(read_file)

    for i in fasta_dict:
        SA, LCP = naive(fasta_dict[i])
        for j in fastq_dict:
            matches_temp = SA_match(SA, LCP, fasta_dict[i], fastq_dict[j])
            if matches_temp:
                for match in matches_temp: # Each iteration makes a SAM row
                    read_name.append(j)
                    reference_name.append(i)
                    match_index.append(match + 1)
                    CIGARS.append(str(len(fastq_dict[j])) + "M")
                    match_string.append(fastq_dict[j])
    output = (read_name, reference_name, match_index, CIGARS, match_string)
    return output

def print_SAM(SAM):
    for i in range(len(SAM[0])):
        sys.stdout.write(SAM[0][i] + "\t" + SAM[1][i] + "\t" + str(SAM[2][i]) + "\t" + SAM[3][i] + "\t" + SAM[4][i] + "\n")

SAM = matches_to_SAM(args.reads, args.genome)

print_SAM(SAM)