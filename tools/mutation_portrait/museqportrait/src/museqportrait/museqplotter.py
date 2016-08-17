'''
Contains classes for generating PDF portrait plots

@author: cnielsen
'''
import logging
from scipy import stats
import numpy as np

import matplotlib as mpl
mpl.use('Agg') # avoids dependency of plotting on an X-server
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('ytick',**{'labelsize':6})
rc('xtick',**{'labelsize':6})
rc('axes', **{'labelsize': 12, 'titlesize': 10})
rc('figure',**{'subplot.hspace': 0.4, 'subplot.bottom': 0.075, 'subplot.top': 0.9})
rc('figure',**{'dpi': 72})

logger = logging.getLogger(__name__)

class Plotter:
    
    def generateSampleReport(self, sample, output, plot_all):
        istring = "generating report " + output
        logger.info(istring)   
        print istring
        
        # set global report layout
        columns = 2
        rows = 4
        if sample.hasTrinucData():
            rows += 1
        if sample.hasCopyNumberData():
            rows += 2
        gs = gridspec.GridSpec(rows, columns)
        fig = plt.figure(figsize=(8.5,11))
        
        # figure header
        fig.text(0.125, 0.97, sample.getName(), weight="bold", size=10,
                 verticalalignment="top")
        fig.text(0.125, 0.95, sample.getStats(), size=10, 
                 verticalalignment="top")

        # add the components
        r = 0;
        logger.info("plotting probability density")
        # plt.subplot(gs[r, 1]);
        plt.subplot(gs[r, 0]);
        self.probDensity(sample)
        plt.subplot(gs[r, 1]);
        self.probDensityCDF(sample)
        
        logger.info("plotting allele ratio versus probability")
        r += 1
        ax = plt.subplot(gs[r, :]);
        self.alleleRatioVersusProb(ax, sample, sample.getThresholds()[0])
        
        r += 1
        x = 0
        for t in sample.getThresholds():
            logger.info("plotting allele ratio versus coverage for threshold " + str(t))
            ax = plt.subplot(gs[r, x]);
            self.alleleRatioVersusCov(ax, sample, t)
            x += 1

        if sample.hasTrinucData():
            r += 1 # this plot on a new row
            logger.info("plotting trinucleotide context")
            ax = plt.subplot(gs[r, :]);
            self.triNucDensityPlot(ax, sample)
        
        if sample.hasCopyNumberData():
            logger.info("plotting copy number data")
            r += 1 # this plot on a new row 
            ax = plt.subplot(gs[r, :]);
            self.copyNumberPlot(ax, sample, plot_all)
            
            logger.info("plotting loh allelic ratio data")
            r += 1
            ax = plt.subplot(gs[r, :]);
            self.allelicRatioPlot(ax, sample, plot_all)
        
        logger.info("plotting genomic distribution")
        r += 1 # this plot on a new row 
        ax = plt.subplot(gs[r, :]);
        self.posDensity(ax, sample)

        # output as pdf
        logger.info("saving plots to pdf")
        p = PdfPages(output)
        p.savefig()
        p.close()
        
    def alleleRatioVersusCov(self, ax, sample, threshold):
        (covs, ratios) = sample.getAlleleRatioAndCov(threshold)
        self.smoothScatter(ax, covs, ratios)
        plt.xlabel('coverage')
        plt.ylabel('allele ratio')
        points = ax.get_position().get_points()
        (x1, y2) = (points[0][0], points[1,1])
        plt.figtext(x1, y2 + 0.005, "P > " + str(threshold), 
                    horizontalalignment="left", verticalalignment="bottom", 
                    fontsize=8, weight="bold")
        
    def alleleRatioVersusProb(self, ax, sample, threshold):
        (probs, ratios) = sample.getAlleleRatioAndProb(threshold)
        self.beanPlot(ax, probs, ratios, 0.1)
        
    def beanPlot(self, ax, xValues, yValues, xBinSize):
        xValues = np.array(xValues); yValues = np.array(yValues)
        xmin = xValues.min(); xmax = xValues.max()
        
        # bin values for bean plotting
        t = xmin + xBinSize
        ySet = []; ySets = []; bins = [xmin]
        for i,v in enumerate(xValues):
            if v < t or v == xmax:
                ySet.append(yValues[i])
            else:
                ySets.append(ySet)
                bins.append(t)
                ySet = [yValues[i]]
                t += xBinSize
        ySets.append(ySet)
      
        w = xBinSize/2 * 0.8
        vmax = []; vvalues = []; xvalues = []
        for i,b in enumerate(bins):
            k = stats.gaussian_kde(ySets[i]) #calculates the kernel density
            s = 1.5*np.std(ySets[i]) # padding to prevent truncation at max and min
            m = k.dataset.min() -s #lower bound of violin
            M = k.dataset.max() +s #upper bound of violin
            x = np.linspace(m, M, 100) # support for violin
            v = k.evaluate(x) #violin profile (density curve)
            # print str(i) + " " + str(b) + " " + str(v.max()) + " " + str(w)
            vmax.append(v.max())
            v = w*v/v.max() #scaling the violin to the available space
            vvalues.append(v)
            xvalues.append(x)
            
            ax.fill_betweenx(x, -v+b, v+b, facecolor='#CAC9CA', edgecolor="#CAC9CA")
            
        # add the data
        hw = w/4
        for i,b in enumerate(bins):
            values = ySets[i]
            k = stats.gaussian_kde(values)
            duplicates = dict()
            for v in values:
                if duplicates.has_key(v): 
                    duplicates[v] += 1
                else:
                    duplicates[v] = 1
            
            # plot the lines
            uniq_values = sorted(duplicates.keys())
            max_count = max(duplicates.values())
            for i,v in enumerate(uniq_values):
                a = max(0.1, float(duplicates[v])/max_count)
                ax.hlines(v, b-hw, b+hw, lw=0.5, color='k', alpha=a)
                
            ax.hlines(np.mean(values), b-w, b+w, lw=1.0, color='k')

        plt.xticks(bins)
        plt.xlabel('probability')
        plt.ylabel('allele ratio')
        plt.ylim(0.0, 1.0)

    def allelicRatioPlot(self, ax, sample, subsample):        
        cTypes = ["HOMD","NLOH","GAIN","DLOH", "ALOH","HET","ASCNA","BCNA","UBCNA"]
        cColors = ["#00FF00","#0000FF","#8B0000","#006400","#006400","#BEBEBE","#FF0000","#BEBEBE","#FF0000"]
        
        points = ax.get_position().get_points()
        (x1,y1) = points[0]
        (x2,y2) = points[1]
        width_in_points = int((x2 - x1) * 8.5 * 72) # 72 dpi, 8.5 in wide, x pos fraction of width
        height_in_points = int((y2 - y1) * 11 * 72)
        yRange = (0.0, 1.0)
        
        print "loh"
        for i,cType in enumerate(cTypes):
            xPositions, logRatios = sample.getAllelicRatiosGenomeWide(subsample, width_in_points, yRange, height_in_points, cType)
            print cColors[i] + " " + cTypes[i] + " " + str(len(logRatios))
            plt.scatter(xPositions, logRatios, s=1, c=cColors[i], marker=',', edgecolor='none')
        # plt.xticks([])
        plt.xlim(0, sum(sample.getChromSizes()))
        plt.ylim(yRange)
        plt.ylabel("allelic ratio")

    def copyNumberPlot(self, ax, sample, plot_all):
        cTypes = ["0","1","2","3","4","5"]
        cColors = ["#00FF00","#006400","#0000FF","#8B0000","#FF0000","#FF0000"]
        
        points = ax.get_position().get_points()
        (x1,y1) = points[0]
        (x2,y2) = points[1]
        width_in_points = int((x2 - x1) * 8.5 * 72) # 72 dpi, 8.5 in wide, x pos fraction of width
        height_in_points = int((y2 - y1) * 11 * 72)
        yRange = (-2.0, 2.0)
        
        print "copy number"
        for i,cType in enumerate(cTypes):
            xPositions, logRatios = sample.getLogRatiosGenomeWide(plot_all, width_in_points, yRange, height_in_points, cType)
            print cColors[i] + " " + cTypes[i] + " " + str(len(logRatios))
            # plt.plot(xPositions, logRatios, 'k.', markersize=1, color=cColors[i])
            plt.scatter(xPositions, logRatios, s=1, c=cColors[i], marker=',', edgecolor='none')
        
        plt.xlim(0, sum(sample.getChromSizes()))
        # plt.xticks([])
        plt.ylim(yRange)
        plt.ylabel("copy number (log ratio)")
        
    def getFilterTypes(self, sample):
        fTypes = []
        if sample.getFilterType() == None:
            # collect all observed filter types in the data
            fTypes = sorted(sample.getProbDensity().keys())
        else:
            # only plot the filtered type
            fTypes = [sample.getFilterType()]
        return fTypes
    
    def probDensity(self, sample):
        pDensity = sample.getProbDensity()
        fTypes = self.getFilterTypes(sample)
        for fType in fTypes:
            pValues = sorted(pDensity[fType].keys())
            counts = [pDensity[fType][p] for p in pValues]
            c = 'k'
            if fType == "PASS":
                c = "#E42927"
            elif fType == "INDL":
                c = "0.5"
            plt.plot(pValues, counts, color = c, linewidth=1.5)
            plt.xlabel('probability')
            plt.ylabel('frequency')
            plt.xlim(sample.getThresholds()[0], 1.0)
            
    def probDensityCDF(self, sample):
        pDensity = sample.getProbDensity()
        fTypes = self.getFilterTypes(sample)
        for fType in fTypes:
            pValues = sorted(pDensity[fType].keys())
            counts = np.array([float(pDensity[fType][p]) for p in pValues])
            cSum = counts.sum()
            counts = counts/cSum
            cdf = np.cumsum(counts)
            c = 'k'
            if fType == "PASS":
                c = "#E42927"
            elif fType == "INDL":
                c = "0.5"
            plt.plot(pValues, cdf, color = c, linewidth=1.5, label=fType)
            plt.xlabel('probability')
            plt.ylabel('CDF')
            plt.xlim(sample.getThresholds()[0], 1.0)
            plt.ylim(0.0, 1.0)
        if len(fTypes) > 1:
            plt.legend(loc = 'lower right', fontsize="xx-small")
        
    def posDensity(self, ax, sample):
        posDensitySets = sample.getGlobalPosDensity()
        numBins = sum([len(l) for l in posDensitySets])
        densities = []; chrLabelPos = []; linePos= []
        
        prevLabelPos = 0.0
        chromNames = sample.getChromNames()
        for i in range(len(chromNames)):
            densities.extend(posDensitySets[i])
            n = len(posDensitySets[i])
            chrLabelPos.append(prevLabelPos + float(n)/2)
            linePos.append(prevLabelPos + n)
            prevLabelPos += n
        
        m = max(densities)
        ax.bar(range(numBins), densities, color='#646464', edgecolor="#646464")
        plt.xlim(0, numBins)
        plt.ylim(0, m)
        plt.ylabel('mutations per Mbp')
        plt.xticks([])
        # plt.xlim(478,670) # chr3
        # plt.xlim(1031,1196)
        
        # chrom labels
        points = ax.get_position().get_points()
        (x1,y1) = points[0]
        (x2,y2) = points[1]
        plotWidth = x2-x1
        chrLabelPos = [x1 + (p/numBins)*plotWidth for p in chrLabelPos]

        for i in range(len(chrLabelPos)):
            plt.figtext(chrLabelPos[i], y1-0.01, chromNames[i], horizontalalignment="center", fontsize=6)
            plt.plot([linePos[i], linePos[i]], [0.0, m], color='#AAAAAA', linewidth=0.5, linestyle="dotted")
        plt.figtext(x1 + plotWidth/2, y1-0.03, "chromosome", horizontalalignment="center", fontsize=10)    
        
        # P-value threshold label
        plt.figtext(x1, y2 + 0.005, "P > " + str(sample.getThresholds()[-1]), 
                    horizontalalignment="left", verticalalignment="bottom", 
                    fontsize=8, weight="bold")     
        
    def pickPlotBoundsOnZScore(self, values, z):
        vMean = np.mean(values)
        vStd = np.std(values)
        maxFromZ = z * vStd + vMean
        minFromZ = -z * vStd + vMean
        lowerBound = max([minFromZ, min(values)])
        upperBound = min([maxFromZ, max(values)])
        return lowerBound, upperBound
        
    def smoothScatter(self, ax, xValues, yValues):
        xValues = np.array(xValues); yValues = np.array(yValues)
        xmin = xValues.min(); xmax = xValues.max()
        ymin = yValues.min(); ymax = yValues.max()
        
        # heatmap
        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        # reformat grid and data values for kernel function
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([xValues, yValues])
        # obtain density estimate
        kernel = stats.gaussian_kde(values)
        # set the covariance_factor, lower means more detail
        kernel.covariance_factor = lambda : .4
        kernel._compute_covariance()
        # reformat density for grid points for plotting
        densities = kernel(positions)
        Z = np.reshape(densities.T, X.shape)
        ax.set_xmargin(0.01)
        ax.set_ymargin(0.01)
        # aspect='auto' is important to enable different scales on x and y axes
        ax.imshow(np.rot90(Z), cmap='Greys', extent=[xmin,xmax,ymin,ymax], aspect='auto')

        # determine a reasonable plotting range (otherwise can be skewed by outliers)
        xLower, xUpper = self.pickPlotBoundsOnZScore(xValues, 3)
        yLower, yUpper = self.pickPlotBoundsOnZScore(yValues, 3)        
        ax.set_xlim([xLower, xUpper])
        ax.set_ylim([yLower, yUpper])
        
    def triNucDensityPlot(self, ax, sample):
        (mutations, contextSets, countSets) = sample.getTriNucData()
        
        # plot each mutation type in a different colour
        cols = ["#1EBEF0", "#070809", "#E42927", "#CAC9CA", "#A1CE64", "#ECC8C5"]
        x = 0
        for i,counts in enumerate(countSets):
            xvalues = [x+j for j in range(len(counts))]
            ax.bar(xvalues, counts, color=cols[i], edgecolor=cols[i])
            x = i*len(counts) + len(counts)
        plt.xlim(0, xvalues[-1])
        plt.ylabel('frequency')
        
        # trinucleotide context labels
        barwidth = 1.0   
        contexts = [item for sublist in contextSets for item in sublist] # flatten the context labels
        ax.set_xticks([x + barwidth/2 for x in range(len(contexts))])
        ax.set_xticklabels(contexts, rotation="vertical")
        
        # mutation labels
        points = ax.get_position().get_points()
        (x1,y1) = points[0]
        (x2,y2) = points[1]
        sectionWidth = (x2-x1)/len(mutations)
        mutLabelPos = [(i*sectionWidth) + sectionWidth/2 + x1 for i in range(len(mutations))]
        for i in range(len(mutations)):
            plt.figtext(mutLabelPos[i], y1-0.035, mutations[i], horizontalalignment="center", fontsize=10)
            
        # P-value threshold label
        plt.figtext(x1, y2 + 0.005, "P > " + str(sample.getThresholds()[-1]), 
                    horizontalalignment="left", verticalalignment="bottom", 
                    fontsize=8, weight="bold")
        
    def testPlot(self, ax, sample):
        
        points = ax.get_position().get_points()
        (x1,y1) = points[0]
        (x2,y2) = points[1]
        
        width_in_points = int((x2 - x1) * 8.5 * 72) # 72 dpi, 8.5 in wide, x pos fraction of width
        height_in_points = int((y2 - y1) * 11 * 72)
        
        # colors = ["#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2", "#5e4fa2"]
        colors = ["#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2"]
        xOffset = 0
        for y in range(height_in_points):
            for i in range(len(colors)):
                xValues = np.array(range(0, width_in_points, len(colors)))
                xValues = xValues + i + xOffset
                plt.scatter(xValues, [y] * len(xValues), s=1, c=colors[i], marker=',', edgecolor='none') # "," = pixel - space-filling as you zoom; s = size in points^2
            xOffset += 1
        plt.xlim(0, width_in_points)
        plt.ylim(0, height_in_points)
            
        
        
