from random import randint
import os

check = None
with open('randomCube.am','r') as in_file, open('newFile.am','w') as out_file:
  for line in in_file:
    if check and line[:2] != '@5':
      out_file.write(str(randint(1,9)))
      out_file.write('\n')
    else:
      out_file.write(line)
    if line[:2] == '@4':
      check = True
os.rename('newFile.am','randomCube.am')
print('The material ids in randomCube.am were successfully randomized.')
