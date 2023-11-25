import importlib.metadata

try:
    __version__ = importlib.metadata.version("pySPM")
except importlib.metadata.PackageNotFoundError:
    __version__ = "unknown"
