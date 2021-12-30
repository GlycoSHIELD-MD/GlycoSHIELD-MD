import os
import setuptools

class CleanCommand(setuptools.Command):
    """Custom clean command to remove unnecessary files."""
    user_options=[]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        cwd=os.path.dirname(os.path.abspath(__file__))
        patterns=['*.pyc', '__pycache__', '*.egg-info']
        for p in patterns:
            os.system("find {} -name \'{}\' -exec rm -rf {{}} \\; -print".format(cwd, p))

setuptools.setup(
    cmdclass={'clean': CleanCommand},
)
