---
name: Bug report
about: Create a report to help us improve

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Snippet of code creating the error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Information:**
 - OS: [e.g. windows]
- Python and pySPM version.

Please run the following and attach the result to your issue
```python
import sys
import pySPM
import numpy as np
import scipy
import matplotlib as mpl
print("Python",sys.version)
print("pySPM",pySPM.__version__)
print("numpy",np.__version__)
print("scipy",scipy.__version__)
print("matplotlib", mpl.__version__)
```

**Additional context**
Add any other context about the problem here.
