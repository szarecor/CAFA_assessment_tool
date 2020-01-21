import os
import pytest
from go_prediction.go_prediction import GOPrediction

@pytest.fixture(scope="module")
def test_data_path():
    ''' Provides a single, consistent absolute path to the test_data directory across environments '''
    root_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return "{}/test/test_data/".format(root_path)

def test_simple_file_parsing(test_data_path):
    ''' Tests the read() method of GOPrediction objects '''
    prediction_path = "ZZZ_1_9606.txt".format(test_data_path)

    test = GOPrediction()
    prediction_handle = open(prediction_path, "r")
    test.read(prediction_handle)

    assert test.author == 'ZZZ'
    assert test.model == 1
    assert test.keywords == ['sequence alignment']
    assert test.taxon == 9606

def test_simple_file_parsing_two(test_data_path):

    obo_path = "{}go_cafa3.obo".format(test_data_path)
    prediction_path = "ZZZ_1_9606.txt".format(test_data_path)

    test = GOPrediction()
    prediction_handle = open(prediction_path, "r")
    test.read_and_split_and_write(obo_path, prediction_handle)

    assert test.author == 'ZZZ'
    assert test.model == 1
    assert test.keywords == ['sequence alignment']
    assert test.taxon == 9606

