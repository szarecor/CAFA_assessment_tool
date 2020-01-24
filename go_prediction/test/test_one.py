import os
import pytest
from go_prediction.go_prediction import GOPrediction, Prediction, PredictionEncoder

@pytest.fixture(scope="module")
def test_data_path():
    ''' Provides a single, consistent absolute path to the test_data directory across environments '''
    root_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return "{}/test/test_data/".format(root_path)

def test_split_by_namespace(test_data_path):
    obo_path = "{}go_cafa3.obo".format(test_data_path)
    prediction_path = "{}ZZZ_1_9606.txt".format(test_data_path)
    parser = GOPrediction()

    '''
    
    
    {'molecular function': [
        Prediction(target='T96060020120', go_term='GO:0008270', confidence=0.8), 
        Prediction(target='T96060020120', go_term='GO:0003700', confidence=0.8), 
        Prediction(target='T96060020120', go_term='GO:0003677', confidence=0.7), 
        Prediction(target='T96060020120', go_term='GO:0003676', confidence=0.29), 
        Prediction(target='T96060020120', go_term='GO:0046872', confidence=0.29), 
        Prediction(target='T96060004121', go_term='GO:0016740', confidence=0.01)
    ], 
    'biological process': [
        Prediction(target='T96060020120', go_term='GO:0006351', confidence=0.8), 
        Prediction(target='T96060020120', go_term='GO:0016032', confidence=0.7), 
        Prediction(target='T96060020120', go_term='GO:0006355', confidence=0.49), 
        Prediction(target='T96060004121', go_term='GO:0008152', confidence=0.01)
    ], 
    'cellular component': [
        Prediction(target='T96060020120', go_term='GO:0005634', confidence=0.7), 
        Prediction(target='T96060020120', go_term='GO:0005730', confidence=0.4), 
        Prediction(target='T96060020120', go_term='GO:0005622', confidence=0.29)
    ]
    }
    '''
    with open(prediction_path, "r") as read_handle:
        res = parser.split_predictions_by_namespace(obo_path=obo_path, prediction_handle=read_handle)
        assert 'molecular function' in res.keys()
        assert 'biological process' in res.keys()
        assert 'cellular component' in res.keys()

        assert len(res.get('molecular function')) == 6
        assert len(res.get('biological process')) == 4
        assert len(res.get('cellular component')) == 3



def test_prediction_file_parser(test_data_path):
    prediction_path = "{}ZZZ_1_9606.txt".format(test_data_path)
    parser = GOPrediction()
    with open(prediction_path, "r") as read_handle:
        parser.parse_prediction_file(read_handle)

        assert 'T96060020120' in parser.data.keys()
        assert len(parser.data.get('T96060020120')) == 11

        T96060020120 = parser.data.get('T96060020120')
        assert T96060020120[0].confidence == .8
        assert T96060020120[1].confidence == .8
        assert T96060020120[-1].confidence == .29

        
        assert 'T96060004121' in parser.data.keys()
        assert len(parser.data.get('T96060004121')) == 2
        '''
        {
        'T96060020120': [
            Prediction(target='T96060020120', go_term='GO:0008270', confidence=0.8), 
            Prediction(target='T96060020120', go_term='GO:0003700', confidence=0.8), 
            Prediction(target='T96060020120', go_term='GO:0006351', confidence=0.8), 
            Prediction(target='T96060020120', go_term='GO:0003677', confidence=0.7), 
            Prediction(target='T96060020120', go_term='GO:0005634', confidence=0.7), 
            Prediction(target='T96060020120', go_term='GO:0016032', confidence=0.7), 
            Prediction(target='T96060020120', go_term='GO:0006355', confidence=0.49), 
            Prediction(target='T96060020120', go_term='GO:0005730', confidence=0.4), 
            Prediction(target='T96060020120', go_term='GO:0005622', confidence=0.29), 
            Prediction(target='T96060020120', go_term='GO:0003676', confidence=0.29), 
            Prediction(target='T96060020120', go_term='GO:0046872', confidence=0.29)
        ], 
            
        'T96060004121': [
            Prediction(target='T96060004121', go_term='GO:0008152', confidence=0.01), 
            Prediction(target='T96060004121', go_term='GO:0016740', confidence=0.01)
        ]
        }
        '''

def test_simple_line_parsing():
    test = GOPrediction()
    test_line = "T96060004121 GO:0008152 0.01"
    parsed_line = test.parse_prediction_line(test_line)
    assert parsed_line.target == 'T96060004121'
    assert parsed_line.go_term == 'GO:0008152'
    assert parsed_line.confidence == 0.01

def test_simple_line_parsing_bad_1():
    test = GOPrediction()
    bad_lines = (
        "AUTHOR ARISTOTLE",
        "MODEL 1",
        "KEYWORDS sequence alignment.",
        "END",
        "T96060004121 GO:0008152 0.95 0.01",
    )

    for bad_line in bad_lines:
        parsed_line = test.parse_prediction_line(bad_line)
        assert parsed_line is None


def test_simple_file_parsing(test_data_path):
    ''' Tests the read() method of GOPrediction objects '''
    prediction_path = "{}ZZZ_1_9606.txt".format(test_data_path)

    test = GOPrediction()
    prediction_handle = open(prediction_path, "r")
    test.read(prediction_handle)

    assert test.author == 'ZZZ'
    assert test.model == 1
    assert test.keywords == ['sequence alignment']
    assert test.taxon == 9606

    print("\n=========================")
    print(test.data)


def test_simple_file_parsing_two(test_data_path):

    obo_path = "{}go_cafa3.obo".format(test_data_path)
    prediction_path = "{}ZZZ_1_9606.txt".format(test_data_path)

    test = GOPrediction()
    prediction_handle = open(prediction_path, "r")
    test.read_and_split_and_write(obo_path, prediction_handle)

    assert test.author == 'ZZZ'
    assert test.model == 1
    assert test.keywords == ['sequence alignment']
    assert test.taxon == 9606

