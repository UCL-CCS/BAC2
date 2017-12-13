
class BaseThermostat:
    pass


class AndersenThermostat(BaseThermostat):
    pass


class VelocityRescale(BaseThermostat):
    def __init__(self, *, temperature, frequency):
        self.temperature = temperature
        self.frequency = frequency

