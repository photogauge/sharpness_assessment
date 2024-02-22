#include "sharpness.h"

#include <pybind11/pybind11.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(_python_api, m) {
    m.doc() = R"pbdoc(
        Python wrapper for `sharpness_assessment`.

        This information will be displayed when using `help()`:
        $ python -c "import sharpness_assessment; help(sharpness_python_api)"
    )pbdoc";

    m.def("assess_sharpness", &assess_sharpness, "Takes an image and assess the sharpness");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
