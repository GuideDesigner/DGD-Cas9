
def return_partner_index(structure):
    """
    Input: structure of a sequence in dot-bracket format
    If structure[i] is '.', return 0,
    If structure[i] is '(', search for the corresponding ')', return the index of the partner.
    If structure[i] is ')', search for the corresponding '(', return the index of the partner.
    """
    partner_index = []
    for i in range(len(structure)):
        if structure[i] == '.':
            partner_index.append(0)
        elif structure[i] == '(':
            left_count = 1
            right_count = 0
            for j in range(i+1, len(structure)):
                if structure[j] == '(':
                    left_count += 1
                elif structure[j] == ')':
                    right_count += 1
                if left_count == right_count:
                    partner_index.append(int(j))
                    break
        elif structure[i] == ')':
            right_count = 1
            left_count = 0
            for j in range(i-1, -1, -1):
                if structure[j] == ')':
                    right_count += 1
                elif structure[j] == '(':
                    left_count += 1
                if left_count == right_count:
                    partner_index.append(int(j))
                    break
    return partner_index

def return_partner_base(parner_index, sequence):
    """
    Input: partner_index, sequence
    Return the base of the partner of each base in the structure.
    """
    partner_base = []
    for i in range(len(parner_index)):
        if parner_index[i] == 0:
            partner_base.append('.')
        else:
            partner_base.append(sequence[parner_index[i]-1])
    return partner_base

def return_dimers_array ( sequence, partner_index ):
    """
    Extract dimers, return a array of dimers
    """
    dimers = []
    for i in range(19):
        if partner_index[i] != 0 and partner_index[i+1] != 0 :
            dimers.append(sequence[i:i+2])
            anti_seq=sequence[partner_index[i+1]]+sequence[partner_index[i]]
            dimers.append(anti_seq)
    return dimers

def adjust_partner_index ( partner_index ):
    """
    adjust partner index
    remove connection which has 1 length
    remain only connection between spacer and scaffold
    """
    for i in range(1,19):
        if (partner_index[i-1]==0 and partner_index[i+1]==0) or partner_index[i]<21:
            partner_index[i]=0
    return partner_index

# make mask array
def return_masked_array(input_array, partner_index, size):
    """
    Input: partner_index
    Return the mask array
    if partner_index[i] == 0, mask_array[i] = ' ' 
    """
    mask_array=[]
    for i in range(0,size):
        if partner_index[i] == 0:
            mask_array.append(' ')
        else :
            mask_array.append(partner_index[i])
    return mask_array

def return_connection_length ( partner_index ):
    """
    figure out each consecutive connection length
    """
    length=[]
    count=0
    cnt=0
    for i in range(20):
        if partner_index[i]!=0:
            count+=1
            cnt+=1
        elif partner_index[i]==0 and count!=0:
            length.append(count)
            count=0
    return(length,cnt)

# convert array to string
def convert_array_to_string ( array ):
    """
    Convert array to string
    """
    return ''.join(array)


# strip the string into array
def strip_string_into_array ( string ):
    """
    Strip the string into array
    """
    string = string.strip()
    return list(string)    

# remove if the size of given element is smaller than given size
def remove_smaller_than_size ( array, size ):
    """
    Remove if the size of given element is smaller than given size
    """
    output_array = []
    for i in range(len(array)):
        if (len(array[i]) >= size):
            output_array.append(array[i])
    return output_array