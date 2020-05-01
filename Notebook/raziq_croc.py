import yt as yt

ds = yt.load("~/Data/rei20c1_a0.1667/rei20c1_a0.1667.art")

def _H1_number_density(field, data):
    return data["H1"] * data["thermal_energy"] * data["density"]
ds.add_field(("gas", "H1_number_density"), function=_raziqpressure, units="1/cm**3")
