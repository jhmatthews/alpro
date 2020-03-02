import alpro.models as models 

class Survival:
	def __init__(self, ModelType):
		self.model = ModelType

		if self.model == "1821":
			self.profile_class = models.RussellProfile()

	def init_model(self):
		self.boxes = models.BoxModel(profile=self.profile_class.profile)

	def get_curve(self, random_seed, L, B0, coherence_func):
		self.boxes.create_box_array(L, B0, random_seed, coherence_func) 
		self.propagate(self.boxes)

	# def propagate():