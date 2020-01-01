
# dev.off()

ht=600; wd=1600
# ht=400; wd=400
# ht=300; wd=300
fn='test.png'
family='mono'
str='SHUD'

png(fn, width=wd, height=ht)
par(mar=c(0,0,0,0))
plot(0,0, type='n',  xlim=c(0, wd), ylim=c(0,ht))
text(wd/2, ht/2, str, cex=wd/35, font=2, family=family)
# grid()
dev.off()