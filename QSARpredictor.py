






class Predictor:

    def __init__(self, input, prRmodels, prout):

        self.input = input
        self.prModels = prRmodels
        self.prout = prout

    def predict(self, model):

        return