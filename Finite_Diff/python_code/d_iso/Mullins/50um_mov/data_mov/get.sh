mkdir pdf jpg
mv *.eps eps
echo eps/*.eps | xargs -n1 pstopdf
mv eps/*.pdf pdf

sips -s format jpeg pdf/*.pdf --out jpg/
