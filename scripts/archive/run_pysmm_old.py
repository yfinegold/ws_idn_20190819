{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "hellooooooooooo\n"
     ]
    }
   ],
   "source": [
    "#set user libraries\n",
    "import ipywidgets as widgets\n",
    "from IPython.display import display, HTML, clear_output, Image\n",
    "from ipywidgets import HBox, Label\n",
    "import subprocess\n",
    "# import ipyleaflet\n",
    "\n",
    "import pysmm,os\n",
    "from pysmm.derive_SM import get_map\n",
    "import sys\n",
    "import ee\n",
    "import itertools\n",
    "from IPython.core.display import HTML \n",
    "\n",
    "\n",
    "ee.Initialize()\n",
    "path = os.path.dirname(pysmm.__file__)\n",
    "user = path.split(\"/\")[2]\n",
    "\n",
    "# year1= sys.argv[1]\n",
    "# month1= sys.argv[2]\n",
    "# day1= sys.argv[3]\n",
    "\n",
    "print('hellooooooooooo')\n",
    "\n",
    "# year1= [2018]\n",
    "# month1= [6]\n",
    "# day1= [1]\n",
    "\n",
    "# for a, b, c in itertools.product( year1, month1, day1):\n",
    "\n",
    "#     print a\n",
    "#     print b\n",
    "#     print c\n",
    "\n",
    "#     get_map(minlon , minlat ,maxlon, maxlat,\n",
    "#             '/data/home/' + user + '/',\n",
    "#             sampling=100,\n",
    "#         year=a, month=b, day=c,\n",
    "#        # year=int(sys.argv[1]), month=int(sys.argv[2]), day=int(sys.argv[3]),\n",
    "#             tracknr=None,\n",
    "#             tempfilter=True,\n",
    "#             mask='Globcover',\n",
    "#             masksnow=False,\n",
    "#             overwrite=True)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "env",
   "language": "python",
   "name": "env"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 2
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython2",
   "version": "2.7.15+"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
