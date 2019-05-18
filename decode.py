import sys
import csv
import random
import numpy as np
import matplotlib.pyplot as plt
import math

def decode(ciphertext, has_breakpoint):

    ##Input:
    ##ciphertext: a ciphered text string
    ##has_breakpoint: a boolean--TRUE means there is a breakpoint, FALSE means not

    csvfile_letter=open("./data/letter_probabilities.csv","r")

    csv_reader_letter = csv.reader(csvfile_letter)
    rows_letter= [row for row in csv_reader_letter]
    letter_prob = np.array(rows_letter)
    letter_prob = np.transpose(letter_prob)
    letter_prob = letter_prob.astype(float)

    csvfile_transition=open("./data/letter_transition_matrix.csv","r")

    csv_reader_transition = csv.reader(csvfile_transition)
    rows_transition= [row for row in csv_reader_transition]
    letter_tran_prob = np.array(rows_transition)
    letter_tran_prob = letter_tran_prob.astype(float)

    cipher = ciphertext
    cipher_size=len(ciphertext)

    dic_cipher = {'a':0, 'b':1, 'c':2, 'd':3, 'e':4, 'f':5, 'g':6, 'h':7, 'i':8, 'j':9,
                  'k':10, 'l':11, 'm':12, 'n':13, 'o':14, 'p':15, 'q':16, 'r':17, 's':18, 't':19
                  , 'u':20, 'v':21, 'w':22, 'x':23, 'y':24, 'z':25, ' ':26, '.':27}
    
    csvfile=open("./data/alphabet.csv","r")

    csv_reader = csv.reader(csvfile)
    rows= [row for row in csv_reader]
    alphabet = np.array(rows)
    alphabet = np.transpose(alphabet)

    alphabet_size=len(alphabet)

    alphabet_index=np.zeros(alphabet_size);
    for i in range(alphabet_size):
        alphabet_index[i]=i
        
    alphabet_map=np.zeros(alphabet_size)
    for i in range(alphabet_size):
        alphabet_map[i]=alphabet_index[i]

    ##If there is a breakpoint:
    if has_breakpoint:
        f_1=np.zeros(alphabet_size)
        for i in range(alphabet_size):
            f_1[i]=alphabet_map[i]

        f_2=np.zeros(alphabet_size)
        for i in range(alphabet_size):
            f_2[i]=alphabet_map[i]

        f_1_final=np.zeros(alphabet_size)
        for i in range(alphabet_size):
            f_1_final[i]=alphabet_map[i]

        f_2_final=np.zeros(alphabet_size)
        for i in range(alphabet_size):
            f_2_final[i]=alphabet_map[i]

        cipher_tmp=np.ones(cipher_size)
        cipher_tmp=cipher_tmp.astype(np.str)
        for i in range(cipher_size):
            cipher_tmp[i]=cipher[i]
        cipher_size_tmp=cipher_size

        pivot=int(cipher_size_tmp/2);
        bp_index=int(cipher_size_tmp/2);

        flag_1=0
        flag_2=0
        flag_first=0

        T=10000
        T_accept=400

        ##Identify the interval where breakpoint exists
        while flag_1==0 or flag_2==0:

            cipher_1=np.ones(bp_index-1)
            cipher_1=cipher_1.astype(np.str)
            for i in range(bp_index-1):
                cipher_1[i]=cipher_tmp[i]

            cipher_2=np.ones(cipher_size_tmp-bp_index)
            cipher_2=cipher_2.astype(np.str)           
            for i in range(cipher_size_tmp-bp_index):
                cipher_2[i]=cipher_tmp[i+(bp_index-1)]

            cipher_size_1=len(cipher_1)
            cipher_size_2=len(cipher_2)
            
            Acceptance_state_record_1=np.zeros(T)
            Acceptance_rate_vector_1=np.zeros(T-T_accept)
            Final_state_1=T-T_accept-1
            
            f_1=np.zeros(alphabet_size)
            for i in range(alphabet_size):
                f_1[i]=alphabet_map[i]
            np.random.shuffle(f_1)

            for t in range(T):
                f_1_tmp=np.random.permutation(f_1)
                f_error=f_1_tmp[0:2]
                f_1_prime=np.zeros(alphabet_size)
                for i in range(alphabet_size):
                    f_1_prime[i]=f_1[i]
                for i in range(alphabet_size):
                    if int(f_1_prime[i])==int(f_error[0]):
                        index_0=i
                    if int(f_1_prime[i])==int(f_error[1]):
                        index_1=i

                tmp=f_1_prime[int(index_0)]
                f_1_prime[int(index_0)]=f_1_prime[int(index_1)]
                f_1_prime[int(index_1)]=tmp

                plain_prime=np.zeros(cipher_size_1)
                plain_prime[0]=f_1_prime[dic_cipher[cipher_1[0]]]
                prod_prime=float(letter_prob[int(plain_prime[0])])
                plain_current=np.zeros(cipher_size_1)
                plain_current[0]=f_1[dic_cipher[cipher_1[0]]]
                prod_current=float(letter_prob[int(plain_current[0])])

                posterior=prod_prime/prod_current
                flag=0

                accept_flag=0

                for i in range(cipher_size_1-2):
                    plain_prime[i+1]=f_1_prime[dic_cipher[cipher_1[i+1]]]
                    prod_prime=float(letter_tran_prob[int(plain_prime[i+1]),int(plain_prime[i])])
                    plain_current[i+1]=f_1[dic_cipher[cipher_1[i+1]]]
                    prod_current=float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])
                    if prod_current>0 and prod_prime>0:
                        posterior*=prod_prime/prod_current
                    elif prod_current==0:
                        flag=1;break
                    else: flag=2;break

                a_prob=0
                if flag==0:
                    a_prob=posterior
                    if a_prob>=1:
                        for i in range(alphabet_size):
                            f_1[i]=f_1_prime[i]
                        accept_flag=1
                    else:
                        a_rand=random.random()
                        if a_rand<a_prob:
                            for i in range(alphabet_size):
                                f_1[i]=f_1_prime[i]
                            accept_flag=1
                if flag==1:
                    for i in range(alphabet_size):
                        f_1[i]=f_1_prime[i]
                    a_prob=1
                    accept_flag=1

                Acceptance_state_record_1[t]=accept_flag

                if flag_1==1:
                    break;

                if t>=T/5 and flag_1==0:
                    Acceptance_rate_vector_1[t-T_accept]=np.sum(Acceptance_state_record_1[t-T_accept:t])/T_accept
                    if Acceptance_rate_vector_1[t-T_accept]<0.3 or t==T:
                        flag_1=1
                        for i in range(alphabet_size):
                            f_1_final[i]=f_1[i]
                    if Acceptance_rate_vector_1[t-T_accept]<0.3 or Acceptance_rate_vector_1[t-T_accept]>0.7:
                        Final_state_1=t-T_accept;break

            Acceptance_state_record_2=np.zeros(T)
            Acceptance_rate_vector_2=np.zeros(T-T_accept)
            Final_state_2=T-T_accept-1

            f_2=np.zeros(alphabet_size)
            for i in range(alphabet_size):
                f_2[i]=alphabet_map[i]
            np.random.shuffle(f_2)

            for t in range(T):
                f_2_tmp=np.random.permutation(f_2)
                f_error=f_2_tmp[0:2]
                f_2_prime=np.zeros(alphabet_size)
                for i in range(alphabet_size):
                    f_2_prime[i]=f_2[i]
                for i in range(alphabet_size):
                    if int(f_2_prime[i])==int(f_error[0]):
                        index_0=i
                    if int(f_2_prime[i])==int(f_error[1]):
                        index_1=i

                tmp=f_2_prime[int(index_0)]
                f_2_prime[int(index_0)]=f_2_prime[int(index_1)]
                f_2_prime[int(index_1)]=tmp

                plain_prime=np.zeros(cipher_size_2)
                plain_prime[0]=f_2_prime[dic_cipher[cipher_2[0]]]
                prod_prime=float(letter_prob[int(plain_prime[0])])
                plain_current=np.zeros(cipher_size_2)
                plain_current[0]=f_2[dic_cipher[cipher_2[0]]]
                prod_current=float(letter_prob[int(plain_current[0])])

                posterior=prod_prime/prod_current
                flag=0

                accept_flag=0

                for i in range(cipher_size_2-2):
                    plain_prime[i+1]=f_2_prime[dic_cipher[cipher_2[i+1]]]
                    prod_prime=float(letter_tran_prob[int(plain_prime[i+1]),int(plain_prime[i])])
                    plain_current[i+1]=f_2[dic_cipher[cipher_2[i+1]]]
                    prod_current=float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])
                    if prod_current>0 and prod_prime>0:
                        posterior*=prod_prime/prod_current
                    elif prod_current==0:
                        flag=1;break
                    else: flag=2;break

                a_prob=0
                if flag==0:
                    a_prob=posterior
                    if a_prob>=1:
                        for i in range(alphabet_size):
                            f_2[i]=f_2_prime[i]
                        accept_flag=1
                    else:
                        a_rand=random.random()
                        if a_rand<a_prob:
                            for i in range(alphabet_size):
                                f_2[i]=f_2_prime[i]
                            accept_flag=1
                if flag==1:
                    for i in range(alphabet_size):
                        f_2[i]=f_2_prime[i]
                    a_prob=1
                    accept_flag=1

                Acceptance_state_record_2[t]=accept_flag

                if flag_2==1:
                    break;

                if t>=T/5 and flag_2==0:
                    Acceptance_rate_vector_2[t-T_accept]=np.sum(Acceptance_state_record_2[t-T_accept:t])/T_accept
                    if Acceptance_rate_vector_2[t-T_accept]<0.3 or t==T:
                        flag_2=1
                        for i in range(alphabet_size):
                            f_2_final[i]=f_2[i]
                    if Acceptance_rate_vector_2[t-T_accept]<0.3 or Acceptance_rate_vector_2[t-T_accept]>0.7:
                        Final_state_2=t-T_accept;break

            if flag_1==0 and flag_2==1:
                cipher_size_tmp=cipher_size_1
                cipher_tmp=np.ones(cipher_size_tmp)
                cipher_tmp=cipher_tmp.astype(np.str)
                for i in range(cipher_size_tmp):
                    cipher_tmp[i]=cipher_1[i]
                pivot=pivot-int(cipher_size_1/2)
                if flag_first==0:
                    flag_first=1
            elif flag_1==1 and flag_2==0:
                cipher_size_tmp=cipher_size_2
                cipher_tmp=np.ones(cipher_size_tmp)
                cipher_tmp=cipher_tmp.astype(np.str)
                for i in range(cipher_size_tmp):
                    cipher_tmp[i]=cipher_2[i]
                pivot_prev=pivot
                pivot=pivot+int(cipher_size_2/2)
                if flag_first==0:
                    flag_first=2
                    
            bp_index=int(cipher_size_tmp/2);

        ##Learn the decipher function by MCMC 
        correct_sign=0
        log_prod_previous_1=0
        
        while correct_sign==0:
            correct_sign=1
            if flag_first==2:
                pivot=pivot_prev
            pivot_std=pivot

            ##Learn the decipher function on the left part of breakpoint
            cipher_1=np.ones(pivot)
            cipher_1=cipher_1.astype(np.str)
            for i in range(pivot-1):
                cipher_1[i]=cipher[i]

            cipher_size_1=len(cipher_1)

            f_cand_size_1=8
            count_log_flag_1=0

            t=0
            
            while t<f_cand_size_1:
                f_1=np.zeros(alphabet_size)
                for i in range(alphabet_size):
                    f_1[i]=alphabet_map[i]
                np.random.shuffle(f_1)

                T_accept=500
                count=0
                total_count=0
                Acceptance_state_record_1=np.zeros(T_accept)

                while total_count<=6000:
                    f_1_tmp=np.random.permutation(f_1)
                    f_error=f_1_tmp[0:2]
                    f_1_prime=np.zeros(alphabet_size)
                    for i in range(alphabet_size):
                        f_1_prime[i]=f_1[i]
                    for i in range(alphabet_size):
                        if int(f_1_prime[i])==int(f_error[0]):
                            index_0=i
                        if int(f_1_prime[i])==int(f_error[1]):
                            index_1=i

                    tmp=f_1_prime[int(index_0)]
                    f_1_prime[int(index_0)]=f_1_prime[int(index_1)]
                    f_1_prime[int(index_1)]=tmp

                    plain_prime=np.zeros(cipher_size_1)
                    plain_prime[0]=f_1_prime[dic_cipher[cipher_1[0]]]
                    prod_prime=float(letter_prob[int(plain_prime[0])])
                    plain_current=np.zeros(cipher_size_1)
                    plain_current[0]=f_1[dic_cipher[cipher_1[0]]]
                    prod_current=float(letter_prob[int(plain_current[0])])

                    posterior=prod_prime/prod_current
                    flag=0

                    accept_flag=0

                    iter_bound_1=min(cipher_size_1-2,5000)

                    for i in range(iter_bound_1):
                        plain_prime[i+1]=f_1_prime[dic_cipher[cipher_1[i+1]]]
                        prod_prime=float(letter_tran_prob[int(plain_prime[i+1]),int(plain_prime[i])])
                        plain_current[i+1]=f_1[dic_cipher[cipher_1[i+1]]]
                        prod_current=float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])
                        if prod_current>0 and prod_prime>0:
                            posterior*=prod_prime/prod_current
                        elif prod_current==0:
                            flag=1;break
                        else: flag=2;break

                    a_prob=0
                    if flag==0:
                        a_prob=posterior
                        if a_prob>=1:
                            for i in range(alphabet_size):
                                f_1[i]=f_1_prime[i]
                            accept_flag=1
                        else:
                            a_rand=random.random()
                            if a_rand<a_prob:
                                for i in range(alphabet_size):
                                    f_1[i]=f_1_prime[i]
                                accept_flag=1
                    if flag==1:
                        for i in range(alphabet_size):
                            f_1[i]=f_1_prime[i]
                        a_prob=1
                        accept_flag=1

                    if count<T_accept-1:
                        count=count+1
                        Acceptance_state_record_1[count-1]=accept_flag
                    else:
                        count=count+1
                        Acceptance_state_record_1[count-1]=accept_flag
                        Acceptance_rate_1=np.sum(Acceptance_state_record_1)/T_accept
                        if Acceptance_rate_1==0:
                            break
                        Acceptance_state_record_1=np.zeros(T_accept)
                        count=0

                    total_count=total_count+1

                plain_current_1=np.zeros(cipher_size_1)
                plain_current_1[0]=f_1[dic_cipher[cipher_1[0]]]
                log_prod_current_1=-math.log(float(letter_prob[int(plain_current_1[0])]))
                for i in range(cipher_size_1-2):
                    plain_current_1[i+1]=f_1[dic_cipher[cipher_1[i+1]]]
                    if float(letter_tran_prob[int(plain_current_1[i+1]),int(plain_current_1[i])])>0:
                        log_trans_tmp_1=-math.log(float(letter_tran_prob[int(plain_current_1[i+1]),int(plain_current_1[i])]))
                        log_prod_current_1=log_prod_current_1+log_trans_tmp_1
                    else:
                        t=t-1
                        log_prod_current_1=0
                        break

                if log_prod_current_1!=0:
                    if count_log_flag_1>0 and log_prod_current_1<log_prod_previous_1:
                        f_1_final=np.zeros(alphabet_size)
                        for i in range(alphabet_size):
                            f_1_final[i]=f_1[i]
                        log_prod_previous_1=log_prod_current_1
                    if count_log_flag_1==0 and log_prod_previous_1==0:
                        f_1_final=np.zeros(alphabet_size)
                        for i in range(alphabet_size):
                            f_1_final[i]=f_1[i]
                        log_prod_previous_1=log_prod_current_1
                    count_log_flag_1+=1

                t=t+1

            ##Precisely find the breakpoint by examining the sudden drop of likelihood function
            index_judge_vector=-np.ones(alphabet_size)
            plain_current=np.zeros(cipher_size)
            plain_current[0]=f_1_final[dic_cipher[cipher[0]]]
            log_prod_current=-math.log(float(letter_prob[int(plain_current[0])]))

            min_prob=1;
            
            for i in range(cipher_size-1):
                plain_current[i+1]=f_1_final[dic_cipher[cipher[i+1]]]
                if index_judge_vector[dic_cipher[cipher[i+1]]]==-1:
                    index_judge_vector[dic_cipher[cipher[i+1]]]=plain_current[i+1]
                if float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])>0:
                    if i>=cipher_size_1-2 and float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])<min_prob:
                        pivot=i
                        break
                    if min_prob>float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])]):
                        min_prob=float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])
                    log_trans_tmp=-math.log(float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])]))
                    log_prod_current=log_prod_current+log_trans_tmp
                else:
                    log_prod_current=0
                    pivot=i
                    break

            ##Learn the decipher function on the right part of breakpoint
            cipher_2=np.ones(cipher_size-2-pivot)
            cipher_2=cipher_2.astype(np.str)           
            for i in range(cipher_size-2-pivot):
                cipher_2[i]=cipher[i+pivot+1]

            cipher_size_2=len(cipher_2)

            f_cand_size_2=8
            count_log_flag_2=0

            log_prod_previous_2=0

            t=0

            while t<f_cand_size_2:
                if pivot<pivot_std and flag_first==2:
                    correct_sign=0
                    pivot=pivot_std
                    break
                if pivot<pivot_std and flag_first==1:
                    correct_sign=0
                    pivot=pivot_std
                    break

                f_2=np.zeros(alphabet_size)
                for i in range(alphabet_size):
                    f_2[i]=alphabet_map[i]
                np.random.shuffle(f_2)

                T_accept=500
                count=0
                total_count=0
                Acceptance_state_record_2=np.zeros(T_accept)

                while total_count<=6000:
                    f_2_tmp=np.random.permutation(f_2)
                    f_error=f_2_tmp[0:2]
                    f_2_prime=np.zeros(alphabet_size)
                    for i in range(alphabet_size):
                        f_2_prime[i]=f_2[i]
                    for i in range(alphabet_size):
                        if int(f_2_prime[i])==int(f_error[0]):
                            index_0=i
                        if int(f_2_prime[i])==int(f_error[1]):
                            index_1=i

                    tmp=f_2_prime[int(index_0)]
                    f_2_prime[int(index_0)]=f_2_prime[int(index_1)]
                    f_2_prime[int(index_1)]=tmp

                    plain_prime=np.zeros(cipher_size_2)
                    plain_prime[0]=f_2_prime[dic_cipher[cipher_2[0]]]
                    prod_prime=float(letter_prob[int(plain_prime[0])])
                    plain_current=np.zeros(cipher_size_2)
                    plain_current[0]=f_2[dic_cipher[cipher_2[0]]]
                    prod_current=float(letter_prob[int(plain_current[0])])

                    posterior=prod_prime/prod_current
                    flag=0

                    accept_flag=0

                    iter_bound_2=min(cipher_size_2-2,5000)

                    for i in range(iter_bound_2):
                        plain_prime[i+1]=f_2_prime[dic_cipher[cipher_2[i+1]]]
                        prod_prime=float(letter_tran_prob[int(plain_prime[i+1]),int(plain_prime[i])])
                        plain_current[i+1]=f_2[dic_cipher[cipher_2[i+1]]]
                        prod_current=float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])
                        if prod_current>0 and prod_prime>0:
                            posterior*=prod_prime/prod_current
                        elif prod_current==0:
                            flag=1;break
                        else: flag=2;break

                    a_prob=0
                    if flag==0:
                        a_prob=posterior
                        if a_prob>=1:
                            for i in range(alphabet_size):
                                f_2[i]=f_2_prime[i]
                            accept_flag=1
                        else:
                            a_rand=random.random()
                            if a_rand<a_prob:
                                for i in range(alphabet_size):
                                    f_2[i]=f_2_prime[i]
                                accept_flag=1
                    if flag==1:
                        for i in range(alphabet_size):
                            f_2[i]=f_2_prime[i]
                        a_prob=1
                        accept_flag=1

                    if count<T_accept-1:
                        count=count+1
                        Acceptance_state_record_2[count-1]=accept_flag
                    else:
                        count=count+1
                        Acceptance_state_record_2[count-1]=accept_flag
                        Acceptance_rate_2=np.sum(Acceptance_state_record_2)/T_accept
                        if Acceptance_rate_2==0:
                            break
                        Acceptance_state_record_2=np.zeros(T_accept)
                        count=0

                    total_count=total_count+1

                plain_current_2=np.zeros(cipher_size_2)
                plain_current_2[0]=f_2[dic_cipher[cipher_2[0]]]
                log_prod_current_2=-math.log(float(letter_prob[int(plain_current_2[0])]))
                for i in range(cipher_size_2-2):
                    plain_current_2[i+1]=f_2[dic_cipher[cipher_2[i+1]]]
                    if float(letter_tran_prob[int(plain_current_2[i+1]),int(plain_current_2[i])])>0:
                        log_trans_tmp_2=-math.log(float(letter_tran_prob[int(plain_current_2[i+1]),int(plain_current_2[i])]))
                        log_prod_current_2=log_prod_current_2+log_trans_tmp_2
                    else:
                        log_prod_current_2=0
                        t=t-1
                        break

                if log_prod_current_2!=0:
                    if count_log_flag_2>0 and log_prod_current_2<log_prod_previous_2:
                        f_2_final=np.zeros(alphabet_size)
                        for i in range(alphabet_size):
                            f_2_final[i]=f_2[i]
                        log_prod_previous_2=log_prod_current_2
                    if count_log_flag_2==0 and log_prod_previous_2==0:
                        f_2_final=np.zeros(alphabet_size)
                        for i in range(alphabet_size):
                            f_2_final[i]=f_2[i]
                        log_prod_previous_2=log_prod_current_2
                    count_log_flag_2+=1

                t=t+1

        ##Decipher the ciphertext by two decipher functions        
        plaintext = np.empty(cipher_size,dtype='str')
        for i in range(cipher_size):
            if i<=pivot:
                plaintext[i]=findkey(dic_cipher,f_1_final[dic_cipher[cipher[i]]])
            else:
                plaintext[i]=findkey(dic_cipher,f_2_final[dic_cipher[cipher[i]]])
        plaintext_string=''.join(plaintext)

    ##If there is not a breakpoint:
    else:
        f_cand_size=8
        count_log_flag=0

        f_final=np.zeros(alphabet_size)
        for i in range(alphabet_size):
            f_final[i]=alphabet_map[i]

        t=0

        ##Run f_cand_size times to avoid being trapped in local extreme points
        while t<f_cand_size:

            f=np.zeros(alphabet_size)
            for i in range(alphabet_size):
                f[i]=alphabet_map[i]
            np.random.shuffle(f)

            T_accept=400
            count=0
            total_count=0
            Acceptance_state_record=np.zeros(T_accept)

            log_prod_previous=0
            while total_count<=6000:
                f_tmp=np.random.permutation(f);
                f_error=f_tmp[0:2]
                f_prime=np.zeros(alphabet_size)
                for i in range(alphabet_size):
                    f_prime[i]=f[i]
                for i in range(alphabet_size):
                    if int(f_prime[i])==int(f_error[0]):
                        index_0=i
                    if int(f_prime[i])==int(f_error[1]):
                        index_1=i

                tmp=f_prime[int(index_0)]
                f_prime[int(index_0)]=f_prime[int(index_1)]
                f_prime[int(index_1)]=tmp
                
                plain_prime=np.zeros(cipher_size)
                plain_prime[0]=f_prime[dic_cipher[cipher[0]]]
                prod_prime=float(letter_prob[int(plain_prime[0])])
                plain_current=np.zeros(cipher_size)
                plain_current[0]=f[dic_cipher[cipher[0]]]
                prod_current=float(letter_prob[int(plain_current[0])])

                posterior=prod_prime/prod_current
                flag=0

                accept_flag=0

                iter_bound=min(cipher_size-2,5000)

                for i in range(iter_bound):
                    plain_prime[i+1]=f_prime[dic_cipher[cipher[i+1]]]
                    prod_prime=float(letter_tran_prob[int(plain_prime[i+1]),int(plain_prime[i])])
                    plain_current[i+1]=f[dic_cipher[cipher[i+1]]]
                    prod_current=float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])
                    if prod_current>0 and prod_prime>0:
                        posterior*=prod_prime/prod_current
                    elif prod_current==0:
                        flag=1;break
                    else: flag=2;break

                a_prob=0
                if flag==0:
                    a_prob=posterior
                    if a_prob>=1:
                        for i in range(alphabet_size):
                            f[i]=f_prime[i]
                        accept_flag=1
                    else:
                        a_rand=random.random()
                        if a_rand<a_prob:
                            for i in range(alphabet_size):
                                f[i]=f_prime[i]
                            accept_flag=1
                if flag==1:
                    for i in range(alphabet_size):
                        f[i]=f_prime[i]
                    a_prob=1
                    accept_flag=1

                if count<T_accept-1:
                    count=count+1
                    Acceptance_state_record[count-1]=accept_flag
                else:
                    count=count+1
                    Acceptance_state_record[count-1]=accept_flag
                    Acceptance_rate=np.sum(Acceptance_state_record)/T_accept
                    if Acceptance_rate==0:
                        break
                    Acceptance_state_record=np.zeros(T_accept)
                    count=0

                total_count=total_count+1

            plain_current=np.zeros(cipher_size)
            plain_current[0]=f[dic_cipher[cipher[0]]]
            log_prod_current=-math.log(float(letter_prob[int(plain_current[0])]))
            for i in range(cipher_size-1):
                plain_current[i+1]=f[dic_cipher[cipher[i+1]]]
                if float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])])>0:
                    log_trans_tmp=-math.log(float(letter_tran_prob[int(plain_current[i+1]),int(plain_current[i])]))
                    log_prod_current=log_prod_current+log_trans_tmp
                else:
                    log_prod_current=0
                    t=t-1
                    break

            if log_prod_current!=0:
                if count_log_flag>0 and log_prod_current<log_prod_previous:
                    f_final=np.zeros(alphabet_size)
                    for i in range(alphabet_size):
                        f_final[i]=f[i]
                    log_prod_previous=log_prod_current
                if count_log_flag==0:
                    f_final=np.zeros(alphabet_size)
                    for i in range(alphabet_size):
                        f_final[i]=f[i]
                    log_prod_previous=log_prod_current
                count_log_flag+=1

            t=t+1
 
        plaintext = np.empty(cipher_size,dtype='str')
        for i in range(cipher_size):
            plaintext[i]=findkey(dic_cipher,f_final[dic_cipher[cipher[i]]])
        plaintext_string=''.join(plaintext)
        
    return plaintext_string

##Find key based on value in dictionary structure
def findkey(dictionary,target_val):
    for key,val in dictionary.items():
        if val == target_val:
            return key
