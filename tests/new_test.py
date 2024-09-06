import sys
import os

print("Python executable:", sys.executable)
print("Current working directory:", os.getcwd())
print("Python path:")
for path in sys.path:
    print(path)

try:
    from plotly_ggtranscript import geom_range, geom_intron, to_intron
    print("Successfully imported plotly_ggtranscript functions")
except ImportError as e:
    print(f"Failed to import plotly_ggtranscript: {e}")

# Rest of your script follows...