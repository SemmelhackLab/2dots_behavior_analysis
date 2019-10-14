import os
import csv
from itertools import izip_longest
import numpy as np

def expt_reader(eyefile):
    '''
    :param eyefile: csv file containing the left and right eye angles
    :return: the left and right eye angles
    '''
    # READ EYE POSITION from the eye tracker

    expt_results = [] # a list to store the angles

    with open(eyefile) as csvfile:

        filename, file_extension = os.path.splitext(eyefile) # split the filename
        readCSV = csv.reader(csvfile, delimiter=',') # declare the csv file separated by comma
        headers = readCSV.next() # store the headers of each column of interest
        fish = headers.index('Fish') # first column with title left
        side = headers.index('Prey Side')
        size = headers.index('Prey Size') # first column with title left
        resp = headers.index('Response time') # first column with title left
        dur = headers.index('Duration') # first column with title left

        for row in readCSV: # read the data of each rows of each column
            expt_results.append([row[fish], row[side], row[size], row[resp], row[dur]])

    return expt_results


def match(fishes, options):
    match_samples = []
    # loop all through samples
    for fish in fishes:
        opp_ID = 0 # ID for finding the opposite side
        same_ID = 0 # ID for finding the same side
        count = -1
        for opt in options:
            count += 1

            # CHECK IF THE SAME FISH
            if fish[0].split('_')[0] != opt[0].split('_')[0]:
                continue
            # Get the fish prey side and size
            side = fish[1]
            size = fish[2]

            # Get the candidate fish prey side and size
            if float(opt[1]) < 0:
                cand_side = "l"
                cand_size = opt[2]
            elif float(opt[1]) > 0:
                cand_side = "r"
                cand_size = opt[2]

            if side != cand_side and opp_ID == 0:
                opp_ID = 1 # matched the condition for opposite side
                match_opp_index = count # save the index for options

            if side == cand_side and same_ID == 0:
                same_ID = 1
                match_same_index = count # save the index for options

            if same_ID == 1 and opp_ID == 1:
                #print fish[0], options[match_opp_index][0], options[match_same_index][0]
                # match_samples.append([fish, options[match_opp_index], options[match_same_index]])
                match_samples.append(fish)
                match_samples.append(options[match_opp_index])
                match_samples.append(options[match_same_index])

                if match_same_index < match_opp_index:
                    options.pop(match_opp_index)
                    options.pop(match_same_index)
                else:
                    options.pop(match_same_index)
                    options.pop(match_opp_index)

                #print "1", options[match_same_index]

                break

    return match_samples


def match2(fishes, options):
    match_samples = []
    # loop all through samples
    for fish in fishes:
        opp_ID = 0 # ID for finding the opposite side
        same_ID = 0 # ID for finding the same side
        count = -1
        for opt in options:
            count += 1

            # CHECK IF THE SAME FISH
            if fish[0].split('_')[0] != opt[0].split('_')[0]:
                continue
            # Get the fish prey side and size
            side = fish[1]
            size = fish[2]

            # Get the candidate fish prey side and size
            if float(opt[1]) < 0:
                cand_side = "l"
                cand_size = opt[2]
            elif float(opt[1]) > 0:
                cand_side = "r"
                cand_size = opt[2]

            if side != cand_side and opp_ID == 0:
                opp_ID = 1 # matched the condition for opposite side
                match_opp_index = count # save the index for options

            if side == cand_side and same_ID == 0:
                same_ID = 1
                match_same_index = count # save the index for options

            if opp_ID == 1:
                # match_samples.append([fish, options[match_opp_index], options[match_same_index]])
                match_samples.append(fish)
                match_samples.append(options[match_opp_index])
                options.pop(match_opp_index)
                break

    return match_samples

dir = "D:\\Semmelhack lab\\002 ANALYSIS\\2p_2dots\\effect of target to competitor\\"
dot = expt_reader(dir + "Summary of 1st PC_1dot_v4.csv")
dots = expt_reader(dir + "Summary of 1st PC_2dots_v4.csv")

A = match2(dot, dots)

results = izip_longest(np.array(A)[:,0][:], np.array(A)[:,1][:], np.array(A)[:,2][:], np.array(A)[:,3][:], np.array(A)[:,4][:])
header = izip_longest(['Fish'], ['Prey Side'], ['Prey Size'], ['Response time'], ['Duration'])

with open(dir + "Summary of all 1st PC_opp" + '.csv', 'wb') as myFile:
     #with open(dir_output + 'Velocity_Acceleration_' + filename + '.csv', 'wb') as myFile:
    wr = csv.writer(myFile, delimiter=',')
    for head in header:
        wr.writerow(head)
    for rows in results:
        wr.writerow(rows)
