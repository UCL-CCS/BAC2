from distutils.version import StrictVersion


class Versioned:
    """Abstract base class for any class that can have a version number.

    """
    @property
    def version(self):
        try:
            return self._version
        except AttributeError:
            return None

    @version.setter
    def version(self, new_version):
        if new_version is None:
            self._version = None
        else:
            self._version = StrictVersion(str(new_version))
