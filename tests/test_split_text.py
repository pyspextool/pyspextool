from pyspextool.utils.split_text import split_text

str = 'This flat was created by scaling the files flat-00017.a.fits, '+ \
  'flat-00018.a.fits, flat-00019.a.fits, flat-00020.a.fits, '+ \
  'flat-00021.a.fits to a common median flux value and then median '+ \
  'combining the scaled images.  The variance is given by '+ \
  '(MAD^2/nimages) where MAD is the median absolute deviation and '+ \
  'nimages is the number of input images.  The zeroth bit of pixels '+ \
  'generated from data with values greater than LINCORMX are set.'

result = split_text(str,length=68)

assert result == ['This flat was created by scaling the files flat-00017.a.fits, ', 'flat-00018.a.fits, flat-00019.a.fits, flat-00020.a.fits, ', 'flat-00021.a.fits to a common median flux value and then median ', 'combining the scaled images.  The variance is given by ', '(MAD^2/nimages) where MAD is the median absolute deviation and ', 'nimages is the number of input images.  The zeroth bit of pixels ', 'generated from data with values greater than LINCORMX are set.']

str = 'This flat was created by combining twenty images.'
result = split_text(str,length=20)

assert result == ['This flat was ', 'created by ', 'combining twenty ',
                  'images.']



