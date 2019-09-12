#import matplotlib libary
import matplotlib.pyplot as plt

#define plot size in inches (width, height) & resolution(DPI)
fig = plt.figure(figsize=(4, 5), dpi=100)

#http://matplotlib.sourceforge.net/api/figure_api.html#matplotlib.figure
ax1 = fig.add_subplot(111)

#define font size
plt.rc("font", size=14)

#define some data
x = [1,2,3,4]
y = [20, 21, 20.5, 20.8]

#error data
y_error = [0.12, 0.13, 0.2, 0.1]

#plot data
ax1.plot(x, y, linestyle="dashed", marker="o", color="green")

#plot only errorbars
ax1.errorbar(x, y, yerr=y_error, linestyle="None", marker="None", color="green")

#configure  X axes
ax1.set_xlim(0.5,4.5)
ax1.set_xticks([1,2,3,4])

#configure  Y axes
ax1.set_ylim(19.8,21.2)
ax1.set_yticks([20, 21, 20.5, 20.8])

#labels
ax1.set_xlabel("this is X",size=10)
ax1.set_ylabel("this is Y", size=10)

#title
plt.title("Simple plot", size=20)

#adjust plot
plt.subplots_adjust(left=0.19)

#show plot
plt.show()
