__author__ = 'Mario'


continzi = [('c1', 'c2'), ('c3', 'c4') ,('c3', 'c5'), ('c1', 'c4'), ('c4', 'c5')]

#formacija puta
def getLast(put):
    return put[-1]

def getFirst(put):
    return put[0]
"""
r = []
put = []
for i in continzi: ##sad su duljine 2
    put.append([i[0], i[1]])
    put.append([i[1], i[0]])

put2 = []
for i in range(len(put)):
    for j in range(len(put)):
        if getLast(put[i]) == getFirst(put[j]) and getLast(put[j]) not in put[i]:
            put2.append(put[i] + [getLast(put[j])])

print(put2)
put3 = []
for i in range(len(put2)):
    for j in range(len(put)):
        if getLast(put2[i]) == getFirst(put[j]) and getLast(put[j]) not in put2[i]:
            put3.append(put2[i] + [getLast(put[j])])
print(put3)

put4 = []
for i in range(len(put3)):
    for j in range(len(put)):
        if getLast(put3[i]) == getFirst(put[j]) and getLast(put[j]) not in put3[i]:
            put4.append(put3[i] + [getLast(put[j])])
print(put4)

put5 = []
for i in range(len(put4)):
    for j in range(len(put)):
        if getLast(put4[i]) == getFirst(put[j]) and getLast(put[j]) not in put4[i]:
            put5.append(put4[i] + [getLast(put[j])])
print(put5)
"""
###idemo ucitat continge
#znaci dobit cu [(ctg1, ctg2), (r1, r0), (r1, ctg0)]

def connect_into_fasta(reads):
    #length = 0
    path0 = []
    paths1 = []
    paths2 = []
    # turn input into a list

    for i in reads:
        path0.append([i[0], i[1]])
        path0.append([i[1], i[0]])

    paths1 = path0[:]
    while True: # maybe while lenght < expected_length
        for i in range(len(paths1)):
            for j in range(len(path0)):
                if getLast(paths1[i]) == getFirst(path0[j]) and getLast(path0[j]) not in paths1[i]:
                    paths2.append(paths1[i] + [getLast(path0[j])])
        if len(paths2) == 0:
            break
        paths1 = paths2
        paths2 = []

    result = []
    for i in paths1:
        if i not in result:
            if i[::-1] not in result:
                result.append(i)
    return result

print(connect_into_fasta(continzi))

