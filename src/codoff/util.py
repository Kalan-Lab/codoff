import os
import traceback

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
        """
        Function from: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
        Call in a loop to create terminal progress bar
        @params:
                iteration   - Required  : current iteration (Int)
                total       - Required  : total iterations (Int)
                prefix      - Optional  : prefix string (Str)
                suffix      - Optional  : suffix string (Str)
                decimals    - Optional  : positive number of decimals in percent complete (Int)
                length      - Optional  : character length of bar (Int)
                fill        - Optional  : bar fill character (Str)
                printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total:        
                print()

def parseCDSCoord(str_gbk_loc):
        try:
                start = None
                end = None
                direction = None
                all_coords = []
                is_multi_part = False
                if not 'join' in str(str_gbk_loc) and not 'order' in str(str_gbk_loc):
                        start = min([int(x.strip('>').strip('<')) for x in
                                                 str(str_gbk_loc)[1:].split(']')[0].split(':')]) + 1
                        end = max([int(x.strip('>').strip('<')) for x in
                                           str(str_gbk_loc)[1:].split(']')[0].split(':')])
                        direction = str(str_gbk_loc).split('(')[1].split(')')[0]
                        all_coords.append([start, end, direction])
                elif 'order' in str(str_gbk_loc):
                        is_multi_part = True
                        all_starts = []
                        all_ends = []
                        all_directions = []
                        for exon_coord in str(str_gbk_loc)[6:-1].split(', '):
                                ec_start = min(
                                        [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                                ec_end = max(
                                        [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                                ec_direction = exon_coord.split('(')[1].split(')')[0]
                                all_starts.append(ec_start)
                                all_ends.append(ec_end)
                                all_directions.append(ec_direction)
                                all_coords.append([ec_start, ec_end, ec_direction])
                        assert (len(set(all_directions)) == 1)
                        start = min(all_starts)
                        end = max(all_ends)
                        direction = all_directions[0]
                else:
                        is_multi_part = True
                        all_starts = []
                        all_ends = []
                        all_directions = []
                        for exon_coord in str(str_gbk_loc)[5:-1].split(', '):
                                ec_start = min(
                                        [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')]) + 1
                                ec_end = max(
                                        [int(x.strip('>').strip('<')) for x in exon_coord[1:].split(']')[0].split(':')])
                                ec_direction = exon_coord.split('(')[1].split(')')[0]
                                all_starts.append(ec_start)
                                all_ends.append(ec_end)
                                all_directions.append(ec_direction)
                                all_coords.append([ec_start, ec_end, ec_direction])
                        assert (len(set(all_directions)) == 1)
                        start = min(all_starts)
                        end = max(all_ends)
                        direction = all_directions[0]
                return(all_coords, start, end, direction, is_multi_part)
        except Exception as e:
                raise RuntimeError(traceback.format_exc())
