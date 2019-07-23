import pandas as pd
import numpy as np

def createVectorScores(line_bins:list, start_zero:bool):
    """
    Create list of values from bins frequencies list. give only a list with the frequencies. The method assumes that the first position correspond the score 0. If start_zero is False, the score zero wasn't considered.

    :param line_bins: list of frequencies
    :param start_zero: to know if the score 0 is considered or not

    :type line_bins: int 
    :type start_zero: int 


    :return: array with the proteins ids
    :rtype: array[int]

    """
    bins_value = 0
    n_repeats = 0
    list_Values = []
    if start_zero == False:
        bins_value = 1
        line_bins.pop(0)

    for n_repeat in line_bins:
        array_bins_aux = np.repeat(bins_value, n_repeat)
        list_Values = np.append(list_Values, array_bins_aux)
        bins_value +=1

    list_Values = list_Values.astype(int)
    return list_Values

def countNumberofScores(list_frequencys:list):
    """
    count the number of scores (not the sum of them but the quantity)

    :param list_frequencys: list of frequencies

    :type line_bins: list(1,4,3,...)


    :return: quantity of elements
    :rtype: int

    """

    list_without_zeros = [element for element in list_frequencys if element != 0]
    qty_elements = np.sum(list_without_zeros)
    return qty_elements


def createVecSizeOfBins(max_score:int, size_bins:int):
    """
    calculate the vector of bins used in case you want a specific size of bins

    :param max_score: highest score that you can have in the bin
    :param size_bins: size of the bins

    :type max_score: int
    :type size_bins: int


    :return: quantity of elements
    :rtype: int

    """
    array_bins_config = np.arange(0, max_score, size_bins)
    if max_score % size_bins != 0:
        array_bins_config = np.concatenate((array_bins_config, [max_score]))
    return array_bins_config

def createVecNumberOfBins(max_score:int, number_of_bins:int):
    """
    calculate the vector of bins used in case you want a specific number of bins

    :param max_score: highest score that you can have in the bin
    :param number_of_bins: number of bins

    :type max_score: int
    :type number_of_bins: int


    :return: vector of bins
    :rtype: array[int]

    """

    size_each_bins = max_score / number_of_bins

    array_bins_config = np.arange(0, max_score, size_each_bins)
    max_score_bin = np.max(array_bins_config)

    if max_score % max_score_bin != 0:
        #array_bins_config.append(max_score)
        array_bins_config = np.concatenate((array_bins_config, [max_score]))
    return array_bins_config

def createHeaderBinsSize(vec_histo_values:np.array):
    """
    create the header for the dataset with the size of bins 

    :param vec_histo_values: vec with the intervals of each bins

    :type vec_histo_values: numPy array

    :return: vec with the labels
    :rtype: int

    """
    vec_designation = ['interaction_ID']
    qty_elements = len(vec_histo_values) - 1
    aux_count = 0
    while aux_count < qty_elements:
        element_start = vec_histo_values[aux_count]
        index_element_end = aux_count + 1
        element_end = vec_histo_values[index_element_end]

        bins_interval_start = float("{0:.3f}".format(element_start))
        bins_interval_end = float("{0:.3f}".format(element_end))
        string_bin_designation = 'bin_' + str(bins_interval_start) + '_' + str(bins_interval_end)
        vec_designation.append(string_bin_designation)

        aux_count += 1

    vec_designation.append('label')
    return vec_designation


def createHeaderBinsNumber(number_of_bins:int):
    """
    create the header for the dataset with the number of bins

    :param number_of_bins: number of bins

    :type number_of_bins: int

    :return: vec with the labels
    :rtype: int

    """
    vec_designation = ['interaction_ID']
    aux = 0
    while aux < number_of_bins:
        string_designation_bins = 'bins_' + str(aux)
        vec_designation.append(string_designation_bins)
        aux += 1
    vec_designation.append('label')
    return vec_designation

def dataTreatmentVecBins(dataframe_csv_base:pd.DataFrame, save_zeros:bool):
    """
    split the data into three vectors
    :note: vec_ids = ids of the interactions
    :note: vec_labels = labels of the interactions
    :note: matrix_frequencies = matric of vectors with the frequencies of each score

    :param number_of_bins: number of bins

    :type number_of_bins: int

    :return: vec with the labels
    :rtype: int

    """
    vec_ids = []
    vec_labels = []
    matrix_frequencies = []

    for index, row in dataframe_csv_base.iterrows():
        id_interaction = row.Interaction_ID
        label_value = row.Label
        list_values = row.values.tolist()
        #Remove the id
        list_values.pop(0)
        #Remove the label
        list_values.pop()
        list_frequencies_hist = createVectorScores(list_values, save_zeros)
        qty_elements = countNumberofScores(list_values)
        qty_scores = len(list_frequencies_hist)
        assert qty_scores == qty_elements

        vec_ids.append(id_interaction)
        vec_labels.append(label_value)
        matrix_frequencies.append(list_frequencies_hist)
        print('N row treated {0}'.format(len(vec_labels)))

    return vec_ids, vec_labels, matrix_frequencies


def constructSizeOfBins(vec_ids:list, vec_labels:list, matrix_frequencies:list, max_score:int, size_of_bins:int):
    """
    This method create the vec of bins score based on the size of bins
    :param vec_ids = ids of the interactions
    :param vec_labels = labels of the interactions
    :param matrix_frequencies = matrix of vectors with the frequencies of each score
    :param max_score = max score found in the dataset based
    :param size_of_bins = size of the bins

    :type vec_ids: list
    :type vec_labels: list
    :type matrix_frequencies: list
    :type max_score: int
    :type size_of_bins: int

    :return: dataframe with the bins
    :rtype: Dataframe

    """
    vec_histo_size = createVecSizeOfBins(max_score, size_of_bins)
    data_bins_treated = []
    for counter, vec_frequencies in enumerate(matrix_frequencies):
        interaction_id = vec_ids[counter]
        label_value = vec_labels[counter]

        histogram_vec_dist = np.histogram(vec_frequencies, bins = vec_histo_size)
        data_histo = histogram_vec_dist[0].flatten()
        data_histo = data_histo.tolist()
        data_histo.append(label_value)
        data_histo.insert(0, interaction_id)
        data_bins_treated.append(data_histo)
    vec_labels = createHeaderBinsSize(histogram_vec_dist[1])
    print(vec_labels)
    dataframe_data = pd.DataFrame(data=data_bins_treated, columns=vec_labels)

    return dataframe_data

def constructNumberOfBins(vec_ids:list, vec_labels:list, matrix_frequencies:list, max_socre:int, number_of_bins:int):
    """
    This method create the vec of bins score based on the number of bins

    :note: the number of bins correspond to the max score/number of bins


    :note vec_ids = ids of the interactions
    :note vec_labels = labels of the interactions
    :note matrix_frequencies = matrix of vectors with the frequencies of each score
    :note max_score = max score found in the dataset based
    :note number_of_bins = size of the bins

    :type vec_ids: list
    :type vec_labels: list
    :type matrix_frequencies: list
    :type max_score: int
    :type number_of_bins: int

    :return: dataframe with the bins
    :rtype: Dataframe

    """
    vec_number_of_bins = createVecNumberOfBins(max_score, number_of_bins)
    data_bins_treated = []
    for counter, vec_frequencies in enumerate(matrix_frequencies):
        interaction_id = vec_ids[counter]
        label_value = vec_labels[counter]

        histogram_vec = np.histogram(vec_frequencies, bins = vec_number_of_bins)
        data_histo = histogram_vec[0].flatten()
        data_histo = data_histo.tolist()
        data_histo.append(label_value)
        data_histo.insert(0, interaction_id)
        data_bins_treated.append(data_histo)
    vec_labels = createHeaderBinsNumber(number_of_bins)
    dataframe_data = pd.DataFrame(data=data_bins_treated, columns=vec_labels)
    print(dataframe_data)
    return dataframe_data

def writeDataframeToCSV(dataframe_bins:pd.DataFrame, file_write:str):
    """
    Save the dataframe to a csv


    :note dataframe_bins = dataframe that contain the bins scores
    :note file_write = name of the file


    :type dataframe_bins: Dataframe
    :type file_write: str


    """
    dataframe_bins.to_csv(index=False, sep=',', path_or_buf=file_write)



path_csv = 'bins_base_573_SB1_ZERO.csv'

dataframe_csv_base = pd.read_csv(filepath_or_buffer=path_csv, delimiter=',')
print(dataframe_csv_base)


dataframe_bins_data = pd.DataFrame()
dataframe_bins_bins_values = pd.DataFrame()

max_score = 573

#take the first 100 rows
dataframe_csv_base = dataframe_csv_base.head(10)
#print(dataframe_csv_base.shape)
use_zeros = False
vec_ids, vec_labels, matrix_frequencies = dataTreatmentVecBins(dataframe_csv_base, use_zeros)
#max_socre = 573
#number_of_bins = 20

#value = constructSizeOfBins(vec_ids, vec_labels, matrix_frequencies, 573, 10)
#value_number_of_bins = constructNumberOfBins(vec_ids, vec_labels, matrix_frequencies, max_socre, number_of_bins)

vec_number_of_bins = [1,5,10,15,30,50]
vec_size_of_bins = [1,5,10,15,20,50]
#performe the number of bins
for number_of_bins in vec_number_of_bins:
    dataframe_results = constructNumberOfBins(vec_ids, vec_labels, matrix_frequencies, max_score, number_of_bins)
    file_name = 'NB_' + str(number_of_bins) + '.csv'
    writeDataframeToCSV(dataframe_results, file_name)

#perform the size of bins
for size_fo_bins in vec_size_of_bins:
    dataframe_results = constructSizeOfBins(vec_ids, vec_labels, matrix_frequencies, max_score, size_fo_bins)
    file_name = 'SB_' + str(size_fo_bins) + '.csv'
    writeDataframeToCSV(dataframe_results, file_name)
