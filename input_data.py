# Plane deformation state. Boolean variable. If False - the plane stress is used. If True - Plane strain is used
PLANE_STRAIN = False

# Soil type. We can choose soil type or set soil characteristics
SOIL_TYPE = 0
# Soil characteristics
Soil_gamma = 0
Soil_E = 0

# Arch types
ARCH_TYPE = 0
# Arch parameters
ARCH_HEIGHT = 5
ARCH_SPAN = 10
# Thickness of the plane scheme
SCHEME_THICKNESS = 1

# Scale for deformed scheme
SCALE_DEF = 50

# friction coefficient
FRICTION_COEFFICIENT = 0.5
# limit of steps for lcp_solve Lemke
LEMKE_LIMIT_STEPS = 800

# Set the value of accuracy for LCP Lemke solver.
# This value means that how many tightening weight (p parameter) we can left to end the solution
ACCURACY_OF_LCP = 10**(-15)

# accuracy of stitching of nodes
ACCURACY_OF_STITCHING = 1e-4

# background color of GUI graphs
BACKGROUND_COLOR = 'w'

# create Excel file in execution directory and write Lemke tables inside
WRITE_EXCEL = True

