import os
import traceback

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
