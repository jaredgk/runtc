import sys



f = open(str(sys.argv[1]),'r')

cnt = int(sys.argv[2])
tag = str(sys.argv[3])

left_f = open(tag+"_left_lengths.txt",'w')
right_f = open(tag+"_right_lengths.txt",'w')

line_count = 0
leftline = ''
rightline = ''
for line in f:
    la = line.strip().split()
    if line_count%cnt == 0:
        if line_count != 0:
            left_f.write(leftline+'\n')
            right_f.write(rightline+'\n')
        leftline = la[1]
        rightline = la[1]
    leftline += ('\t'+la[3])
    rightline += ('\t'+la[2])
    line_count += 1

left_f.write(leftline+'\n')
right_f.write(rightline+'\n')

