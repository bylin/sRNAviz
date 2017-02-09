from tRNA_position import *
import re
import sys

def position_base(positions, seq):
  for position_index, position in enumerate(positions):
    if position.paired:
      index1, index2 = position.position.split(':')
      index1, index2 = int(index1), int(index2)
      base_pair = "{}:{}".format(seq[index1 - 1].upper(), seq[index2 - 1].upper())
      yield position.sprinzl, base_pair
    else:
      index = int(position.position)
      base = seq[index - 1].upper()
      yield position.sprinzl, base


# read in input
seq = sys.stdin.readline()
ss = sys.stdin.readline()

# generate identities from alignment
positions = annotate_positions(ss)
identities = {sprinzl:base for sprinzl, base in position_base(positions, seq)}

# create single base identities from paired identities
cols = list(filter(lambda x: ':' in x, identities.keys()))
for col in cols:
  pos1, pos2 = col.split(':')
  identities[pos1] = identities[col].split(':')[0]
  identities[pos2] = identities[col].split(':')[1]

# create 3d base pairs
identities['8:14'] = '{}:{}'.format(identities['8'], identities['14'])
identities['9:23'] = '{}:{}'.format(identities['9'], identities['23'])
identities['10:45'] = '{}:{}'.format(identities['10'], identities['45'])
identities['22:46'] = '{}:{}'.format(identities['22'], identities['46'])
identities['15:48'] = '{}:{}'.format(identities['15'], identities['48'])
identities['18:55'] = '{}:{}'.format(identities['18'], identities['55'])
identities['19:56'] = '{}:{}'.format(identities['19'], identities['56'])
identities['26:44'] = '{}:{}'.format(identities['26'], identities['44'])
identities['54:58'] = '{}:{}'.format(identities['54'], identities['58'])

for sprinzl in sorted(identities.keys()): sys.stdout.write('{}\t{}\n'.format(sprinzl, identities[sprinzl]))