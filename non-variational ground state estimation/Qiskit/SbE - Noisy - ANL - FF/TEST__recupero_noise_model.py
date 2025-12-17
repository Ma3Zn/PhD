from qiskit_ibm_runtime.fake_provider import FakeBrisbane
from qiskit_aer.noise import NoiseModel


backend = FakeBrisbane()
noise_model = NoiseModel.from_backend(backend)

print(noise_model)