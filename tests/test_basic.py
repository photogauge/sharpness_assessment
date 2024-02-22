from sharpness_python_api import assess_sharpness

def test_add():
    assess_sharpness("./images/fuzzy.jpg", "./images/output/", 0.5)
