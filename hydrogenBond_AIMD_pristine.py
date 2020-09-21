from math import sqrt, ceil
import math
import csv
from random import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def load_positions(filename,newline=''):
    with open(filename) as csvfile:
        output_data = csv.reader(csvfile, delimiter='\t')
        output = []
        for i in output_data:
            row = i
            output.append(row)
        return output
    
def calculate_distance(p1,p2):
    r_vect = p1 - p2
    r_mag = np.linalg.norm(r_vect)
    return r_mag



NUM_OF_LINES=162 #how many atoms per snapshot

filename = 'IPA.txt' #this file contains both framewrok and IPA positions

#split file into separate snapshots
with open(filename) as fin:
    a=2
    fout = open("output1.txt","w")
    for i,line in enumerate(fin):
      fout.write(line)
      if (i+1)%NUM_OF_LINES == 0:
        fout.close()
        fout = open("output%d.txt"%a,"w")
        a=a+1

    fout.close()



IPA_addedOH=[]
IPA_mu3OH=[]
IPA_both=[]
sum_addedOH=0
sum_mu3OH=0
sum_both=0

# loop all the snapshot files
for snapshot in range (1,a-1):
    file = 'output%d.txt'%snapshot


    output=load_positions(file)



    for i in range(1,len(output)):
        ax=float(output[i][1])
        ay=float(output[i][2])
        az=float(output[i][3])
        output[i][1]=ax
        output[i][2]=ay
        output[i][3]=az

# indentify IPA_O IPA_H mu3O mu3H

    output[127-1]=['IPA_O', output[127-1][1], output[127-1][2], output[127-1][3]]
    output[115-1]=['IPA_O', output[115-1][1], output[115-1][2], output[115-1][3]]
    output[151-1]=['IPA_O', output[151-1][1], output[151-1][2], output[151-1][3]]
    output[139-1]=['IPA_O', output[139-1][1], output[139-1][2], output[139-1][3]]


    output[153-1]=['IPA_H', output[153-1][1], output[153-1][2], output[153-1][3]]
    output[117-1]=['IPA_H', output[117-1][1], output[117-1][2], output[117-1][3]]
    output[129-1]=['IPA_H', output[129-1][1], output[129-1][2], output[129-1][3]]
    output[141-1]=['IPA_H', output[141-1][1], output[141-1][2], output[141-1][3]]


    output[101-1]=['mu3O', output[101-1][1], output[101-1][2], output[101-1][3]]
    output[102-1]=['mu3O', output[102-1][1], output[102-1][2], output[102-1][3]]
    output[103-1]=['mu3O', output[103-1][1], output[103-1][2], output[103-1][3]]
    output[104-1]=['mu3O', output[104-1][1], output[104-1][2], output[104-1][3]]

    output[73-1]=['mu3H', output[73-1][1], output[73-1][2], output[73-1][3]]
    output[74-1]=['mu3H', output[74-1][1], output[74-1][2], output[74-1][3]]
    output[75-1]=['mu3H', output[75-1][1], output[75-1][2], output[75-1][3]]
    output[76-1]=['mu3H', output[76-1][1], output[76-1][2], output[76-1][3]]



    a = len(output)

#conert cartesian coordinate to fractional coordinate
    def get_cartesian_to_fractional_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):

        if angle_in_degrees:
            alpha = np.deg2rad(alpha)
            beta = np.deg2rad(beta)
            gamma = np.deg2rad(gamma)
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        r = np.zeros((3, 3))
        r[0, 0] = 1.0 / a
        r[0, 1] = -cosg / (a * sing)
        r[0, 2] = (cosa * cosg - cosb) / (a * volume * sing)
        r[1, 1] = 1.0 / (b * sing)
        r[1, 2] = (cosb * cosg - cosa) / (b * volume * sing)
        r[2, 2] = sing / (c * volume)
        return r


    r=get_cartesian_to_fractional_matrix(14.834, 14.834, 14.834, 60, 60, 60,angle_in_degrees=True)


    for i in range (1,len(output)):
        output[i][1:4] = np.dot(r, output[i][1:4])
        
        
#group mu3H mu3O IPA_O IPA_H separately
    positions_mu3H = []
    for i in range (1,a):
        if 'mu3H' in output[i]:
            positions_mu3H.append(output[i][1:4])

    positions_mu3O = []
    for i in range (1,a):
        if 'mu3O' in output[i]:
            positions_mu3O.append(output[i][1:4])

    positions_IPA_O = []
    for i in range (1,a):
        if 'IPA_O' in output[i]:
            positions_IPA_O.append(output[i][1:4])

    positions_IPA_H = []
    for i in range (1,a):
        if 'IPA_H' in output[i]:
            positions_IPA_H.append(output[i][1:4])


#duplicate mu3H and mu3O to its periodic position if it is within 3 A of edges

    ratio = 3/14.834


    for i in range (0,len(positions_mu3H)):
        if positions_mu3H[i][0]>1-ratio:
            array= [positions_mu3H[i][0]-1, positions_mu3H[i][1], positions_mu3H[i][2]]
            positions_mu3H.append(array)
        if positions_mu3H[i][0]<ratio:
            array= [positions_mu3H[i][0]+1, positions_mu3H[i][1], positions_mu3H[i][2]]
            positions_mu3H.append(array)
        if positions_mu3H[i][1]>1-ratio:
            array= [positions_mu3H[i][0], positions_mu3H[i][1]-1, positions_mu3H[i][2]]
            positions_mu3H.append(array)
        if positions_mu3H[i][1]<ratio:
            array= [positions_mu3H[i][0], positions_mu3H[i][1]+1, positions_mu3H[i][2]]
            positions_mu3H.append(array)
        if positions_mu3H[i][2]>1-ratio:
            array= [positions_mu3H[i][0], positions_mu3H[i][1], positions_mu3H[i][2]-1]
            positions_mu3H.append(array)
        if positions_mu3H[i][2]<ratio:
            array= [positions_mu3H[i][0], positions_mu3H[i][1], positions_mu3H[i][2]+1]
            positions_mu3H.append(array)

    for i in range (0,len(positions_mu3O)):
        if positions_mu3O[i][0]>1-ratio:
            array= [positions_mu3O[i][0]-1, positions_mu3O[i][1], positions_mu3O[i][2]]
            positions_mu3O.append(array)
        if positions_mu3O[i][0]<ratio:
            array= [positions_mu3O[i][0]+1, positions_mu3O[i][1], positions_mu3O[i][2]]
            positions_mu3O.append(array)
        if positions_mu3O[i][1]>1-ratio:
            array= [positions_mu3O[i][0], positions_mu3O[i][1]-1, positions_mu3O[i][2]]
            positions_mu3O.append(array)
        if positions_mu3O[i][1]<ratio:
            array= [positions_mu3O[i][0], positions_mu3O[i][1]+1, positions_mu3O[i][2]]
            positions_mu3O.append(array)
        if positions_mu3O[i][2]>1-ratio:
            array= [positions_mu3O[i][0], positions_mu3O[i][1], positions_mu3O[i][2]-1]
            positions_mu3O.append(array)
        if positions_mu3O[i][2]<ratio:
            array= [positions_mu3O[i][0], positions_mu3O[i][1], positions_mu3O[i][2]+1]
            positions_mu3O.append(array)



#convert fractional to cartesian coordinate

    def get_fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma, angle_in_degrees=True):
        if angle_in_degrees:
            alpha = np.deg2rad(alpha)
            beta = np.deg2rad(beta)
            gamma = np.deg2rad(gamma)
        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)
        r = np.zeros((3, 3))
        r[0, 0] = a
        r[0, 1] = b * cosg
        r[0, 2] = c * cosb
        r[1, 1] = b * sing
        r[1, 2] = c * (cosa - cosb * cosg) / sing
        r[2, 2] = c * volume / sing
        return r



    r=get_fractional_to_cartesian_matrix(14.834, 14.834, 14.834, 60, 60, 60,angle_in_degrees=True)


    for i in range (0,len(positions_mu3O)):
        positions_mu3O[i] = np.dot(r, positions_mu3O[i])

    for i in range (0,len(positions_mu3H)):
        positions_mu3H[i] = np.dot(r, positions_mu3H[i])

    for i in range (0,len(positions_IPA_O)):
        positions_IPA_O[i] = np.dot(r, positions_IPA_O[i])

    for i in range (0,len(positions_IPA_H)):
        positions_IPA_H[i] = np.dot(r, positions_IPA_H[i])



    total_IPA = len(positions_IPA_H)

    print('total number of IPA: '+str(total_IPA))

    sigma_mix = 2 # cutoff distance of hydrogen bond
    count_addedOH = 0
    count_u3OH = 0
    count_addedOH_u3OH =0
    
    #check hydrogen bond for all the IPA

    for j in range (0, len(positions_IPA_H)):


        distance_u = 0
        distance = 0

        list_u = []
        p3 = positions_IPA_H[j]
        p4 = positions_IPA_O[j]


        for i in range (0, len(positions_mu3O)):
            p1 = positions_mu3O[i]
            distance=calculate_distance(p1,p3)
            if distance<sigma_mix:
                list_u.append(distance)
        for i in range (0, len(positions_mu3H)):
            p2 = positions_mu3H[i]
            distance=calculate_distance(p2,p4)
            if distance<sigma_mix:
                list_u.append(distance)
        if list_u:
            distance_u = float(min(list_u))
        else: 
            distance_u = 0        


        if distance_u != 0: 
            count_u3OH +=1






    print("IPA bonded on mu3OH "+str(count_u3OH))
    total = count_u3OH
    print("total number of hydrogen bond  "+str(total))


    percentage_u3OH = count_u3OH/total_IPA


    print("hrdrogen bond on added u3OH/total number of IPA "+str(percentage_u3OH))



    print("number of snapshot"+str(snapshot))
    print("IPA bonded on added OH "+str(count_addedOH))
    print("IPA bonded on mu3OH "+str(count_u3OH))
    print("IPA bonded to both added OH and mu3OH "+str(count_addedOH_u3OH))
    total = count_u3OH+count_addedOH-count_addedOH_u3OH
#     print("total number of hydrogen bond  "+str(total))
    if total_IPA != 0:
        percentage_addedOH = count_addedOH/total_IPA
        percentage_u3OH = count_u3OH/total_IPA
        percentage_hydrogenbond=total/total_IPA
        IPA_addedOH.append(percentage_addedOH)
        IPA_mu3OH.append(percentage_u3OH)
        IPA_both.append(percentage_hydrogenbond)
        sum_addedOH+=percentage_addedOH
        sum_mu3OH+=percentage_u3OH
        sum_both+=percentage_hydrogenbond
    else:
        a=a-1
        print(a)
        
    print("hrdrogen bond on added OH/total number of IPA "+str(percentage_addedOH))
    print("hrdrogen bond on added u3OH/total number of IPA "+str(percentage_u3OH))
    print("hrdrogen bond /total number of IPA "+str(percentage_hydrogenbond))

#plot
x = []
for i in range(1,117):
    x.append(i)


label2="ratio of IPA bonded to mu3OH:"+str(round(sum_mu3OH/(a-1),2))


plt.plot(x, IPA_mu3OH, 'b', label=label2)

plt.legend(loc="best")
plt.ylabel('ratio')
plt.xlabel('snapshot number')
plt.title("4 IPA in pristine UiO-66 AIMD cutoff:2A")
plt.show()