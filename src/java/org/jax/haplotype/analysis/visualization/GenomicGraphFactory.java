/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.analysis.visualization;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Paint;
import java.awt.geom.Rectangle2D;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.RealValuedBasePairInterval;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.SymbolAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.CombinedRangeXYPlot;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYBarPainter;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.chart.renderer.xy.XYBlockRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.DefaultXYZDataset;
import org.jfree.data.xy.XIntervalSeries;
import org.jfree.data.xy.XIntervalSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.RectangleInsets;

/**
 * My little factory for creating graphs from SNP intervals
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenomicGraphFactory
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            GenomicGraphFactory.class.getName());
    
    /**
     * Creates a plot which contains multiple chromosome histograms
     * @param chromosomeHistograms
     *          the chromosome histograms
     * @param xAxisLabel
     *          the X axis label
     * @param yAxisLabel
     *          the Y axis label
     * @return
     *          the object for the multiple chromosome plots
     */
    public JFreeChart createMultiChromosomeHistogram(
            final List<? extends ChromosomeHistogramValues> chromosomeHistograms,
            final String xAxisLabel,
            final String yAxisLabel)
    {
        // create the y-axis that is shared by all of the chromosomes
        NumberAxis yAxis = new NumberAxis();
        if(yAxisLabel != null)
        {
            yAxis.setLabel(yAxisLabel);
        }
        
        // the combined plot makes all of the chromosomes to share the same
        // y-axis but to have thier own x-axis
        CombinedRangeXYPlot combinedChromosomePlot = new CombinedRangeXYPlot(yAxis);
        combinedChromosomePlot.setGap(4.0);
        
        // iterate through the chromosomes adding a new plot for each one
        for(ChromosomeHistogramValues currChromHist: chromosomeHistograms)
        {
            // use a weight to ensure that the amount of space that the
            // subgraph gets on the x-axis is proportional to the chromosome's
            // extent
            final int currWeight =
                (int)(currChromHist.getExtentInBasePairs() / 1000 + 1);
            
            final NumberAxis currXAxis = new NumberAxis();
            currXAxis.setAutoRangeIncludesZero(false);
            currXAxis.setTickLabelsVisible(false);
            currXAxis.setTickMarksVisible(false);
            currXAxis.setLabel(Integer.toString(
                    currChromHist.getChromosomeNumber()));
            
            final XYPlot currPlot = this.createSnpIntervalHistogramPlot(
                    currChromHist.getIntervals(),
                    currChromHist.getVisualInterval(),
                    currXAxis,
                    null);
            
            combinedChromosomePlot.add(currPlot, currWeight);
        }
        
        JFreeChart multiChromosomeChart = new JFreeChart(combinedChromosomePlot);
        multiChromosomeChart.removeLegend();
        
        return multiChromosomeChart;
    }
    
    /**
     * Create a snp interval histogram without any axes
     * @param intervals
     *          the intervals to use
     * @param startInBasePairs
     *          where should we start the graph?
     * @param extentInBasePairs
     *          how far should the graph extend
     * @param visualInterval
     *          the visual interval to use
     * @param xAxisLabel
     *          the x axis title to use
     * @param yAxisLabel
     *          the y axis title to use
     * @return
     *          the histogram
     */
    public JFreeChart createSnpIntervalHistogram(
            final List<? extends RealValuedBasePairInterval> intervals,
            final long startInBasePairs,
            final long extentInBasePairs,
            final HighlightedSnpInterval visualInterval,
            final String xAxisLabel,
            final String yAxisLabel)
    {
        // create the axes
        NumberAxis xAxis = new NumberAxis();
        xAxis.setAutoRangeIncludesZero(false);
        xAxis.setRange(new Range(
                startInBasePairs,
                startInBasePairs + extentInBasePairs));
        if(xAxisLabel != null)
        {
            xAxis.setLabel(xAxisLabel);
        }
        
        NumberAxis yAxis = new NumberAxis();
        if(yAxisLabel != null)
        {
            yAxis.setLabel(yAxisLabel);
        }
        
        // create the plot
        XYPlot plot = this.createSnpIntervalHistogramPlot(
                intervals,
                visualInterval,
                xAxis,
                yAxis);
        
        // create the final chart
        JFreeChart histogram = new JFreeChart(plot);
        histogram.removeLegend();
        
        return histogram;
    }
    
    /**
     * Create a snp interval histogram without any axes
     * @param intervals
     *          the intervals to use
     * @param startInBasePairs
     *          where should we start the graph?
     * @param extentInBasePairs
     *          how far should the graph extend
     * @param visualInterval
     *          the visual interval to use
     * @return
     *          the histogram
     */
    public JFreeChart createSnpIntervalHistogram(
            final List<? extends RealValuedBasePairInterval> intervals,
            final long startInBasePairs,
            final long extentInBasePairs,
            final HighlightedSnpInterval visualInterval)
    {
        // create the axes
        NumberAxis xAxis = new NumberAxis();
        xAxis.setAutoRangeIncludesZero(false);
        xAxis.setRange(new Range(
                startInBasePairs,
                startInBasePairs + extentInBasePairs));
        NumberAxis yAxis = new NumberAxis();
        
        // hide the axes
        xAxis.setVisible(false);
        yAxis.setVisible(false);
        
        // create the plot
        XYPlot plot = this.createSnpIntervalHistogramPlot(
                intervals,
                visualInterval,
                xAxis,
                yAxis);
        
        // more hiding
        plot.setInsets(new RectangleInsets(0.0, 0.0, 0.0, 0.0));
        
        // create the final chart
        JFreeChart histogram = new JFreeChart(plot);
        histogram.removeLegend();
        
        return histogram;
    }
    
    private XYPlot createSnpIntervalHistogramPlot(
            final List<? extends RealValuedBasePairInterval> intervals,
            final HighlightedSnpInterval visualInterval,
            ValueAxis domainAxis,
            ValueAxis rangeAxis)
    {
        XYDataset dataset = this.createSnpIntervalHistogramData(
                intervals,
                visualInterval);
        XYBarRenderer renderer = new XYBarRenderer()
        {
            /**
             * every serializable is supposed to have one of these
             */
            private static final long serialVersionUID = -7907956889918704042L;

            /**
             * {@inheritDoc}
             */
            @Override
            public Paint getItemPaint(int row, int column)
            {
                // Always render in black
                return Color.BLACK;
            }
        };
        renderer.setShadowVisible(false);
        renderer.setDrawBarOutline(false);
//        renderer.setGradientPaintTransformer(null);
        renderer.setBarPainter(new StandardXYBarPainter());
        renderer.setMargin(0.0);
        return new XYPlot(
                dataset,
                domainAxis,
                rangeAxis,
                renderer);
    }
    
    /**
     * Convert the intervals & values to a dataset that JFreeChart can use
     * @param snpIntervals
     *          the intervals
     * @param visualInterval
     *          the visual interval that we should use
     * @return
     *          the data set
     */
    private XYDataset createSnpIntervalHistogramData(
            final List<? extends RealValuedBasePairInterval> snpIntervals,
            final HighlightedSnpInterval visualInterval)
    {
        XIntervalSeriesCollection dataset = new XIntervalSeriesCollection();
        
        XIntervalSeries series = new XIntervalSeries(
                "Interval Values",
                false,
                true);
        XIntervalSeries highlightSeries = new XIntervalSeries(
                "Highlighted Interval Values",
                false,
                true);
        
        int endIndex = visualInterval.getEndIndex();
        int[] highlightIndices = visualInterval.getIndicesToHighlight();
        for(int i = visualInterval.getStartIndex(); i <= endIndex; i++)
        {
            RealValuedBasePairInterval currSnpInterval = snpIntervals.get(i);
            
            boolean highlight = false;
            for(int j = 0; j < highlightIndices.length; j++)
            {
                if(i == highlightIndices[j])
                {
                    highlight = true;
                    break;
                }
            }
            
            if(highlight)
            {
                highlightSeries.add(
                    currSnpInterval.getStartInBasePairs(),
                    currSnpInterval.getStartInBasePairs(),
                    currSnpInterval.getEndInBasePairs(),
                    currSnpInterval.getRealValue());
            }
            else
            {
                series.add(
                        currSnpInterval.getStartInBasePairs(),
                        currSnpInterval.getStartInBasePairs(),
                        currSnpInterval.getEndInBasePairs(),
                        currSnpInterval.getRealValue());
            }
        }
        
        dataset.addSeries(series);
        dataset.addSeries(highlightSeries);
        return dataset;
    }
    
    /**
     * Create a SNP block graph for the given parameters. This graph
     * will show where the intervals do and don't exist. Interval lists
     * are organized with the X axis label matching the list's key
     * @param snpIntervals
     *          the blocks
     * @param startInBasePairs
     *          the x location to start the graph at
     * @param extentInBasePairs
     *          the extent to use for the graph
     * @param maximumImageBlockCount
     *          the max # of separate image blocks to use
     * @param renderAxes
     *          if true render the axes... otherwise dont
     * @param legendText
     *          the text to use for the legend
     * @param yAxisText
     *          the text to use for the y axis
     * @param trueColor
     *          the color to use for true
     * @param falseColor
     *          the color to use for false 
     * @return
     *          the graph
     */
    public JFreeChart createSnpIntervalGraph(
            final Map<String, ? extends List<? extends BasePairInterval>> snpIntervals,
            final long startInBasePairs,
            final long extentInBasePairs,
            final int maximumImageBlockCount,
            final boolean renderAxes,
            final String legendText,
            final String yAxisText,
            final Color trueColor,
            final Color falseColor)
    {
        XYDataset dataset = snpIntervalsToDataset(
                snpIntervals,
                startInBasePairs,
                extentInBasePairs,
                maximumImageBlockCount);
        
        NumberAxis xAxis = new NumberAxis("SNP Position (Base Pairs)");
        xAxis.setAutoRangeIncludesZero(false);
        xAxis.setRange(new Range(
                startInBasePairs,
                startInBasePairs + extentInBasePairs));
        String[] sortedStrainNames = extractSortedStrainNames(snpIntervals);
        for(int strainIndex = 0; strainIndex < sortedStrainNames.length; strainIndex++)
        {
            LOG.info("Strain Name: " + sortedStrainNames[strainIndex]);
        }
        SymbolAxis yAxis = new SymbolAxis(
                yAxisText == null ? "" : yAxisText,
                sortedStrainNames);
        
        XYPlot plot = new XYPlot(dataset, xAxis, yAxis, null);
        LegendItemCollection items = new LegendItemCollection();
        if(legendText != null)
        {
            items.add(new LegendItem(
                    legendText,
                    null,
                    null,
                    null,
                    new Rectangle2D.Double(-3.0, -3.0, 6.0, 6.0),
                    trueColor,
                    new BasicStroke(),
                    Color.BLACK));
        }
        plot.setFixedLegendItems(items);
        
        XYBlockRenderer r = new XYBlockRenderer();
        SmoothPaintScale ps = new SmoothPaintScale(
                0.0,
                1.0,
                falseColor,
                trueColor);

        r.setPaintScale(ps);
        r.setBlockHeight(1.0);
        r.setBlockWidth(1.0);
        plot.setRenderer(r);
        
        final JFreeChart chart;
        if(renderAxes)
        {
            chart = new JFreeChart(
                    "Identical By State Blocks",
                    JFreeChart.DEFAULT_TITLE_FONT,
                    plot,
                    true);
        }
        else
        {
            xAxis.setVisible(false);
            yAxis.setVisible(false);
            plot.setInsets(new RectangleInsets(0.0, 0.0, 0.0, 0.0));
            chart = new JFreeChart(plot);
        }
        
        chart.setBackgroundPaint(Color.WHITE);
        
        return chart;
    }
    
    /**
     * Pull out the keys as a sorted array
     * @param snpIntervals
     *          the interval map (we only care about the keys)
     * @return
     *          the sorted array
     */
    private String[] extractSortedStrainNames(
            final Map<String, ? extends List<? extends BasePairInterval>> snpIntervals)
    {
        String[] sortedStrainNames = snpIntervals.keySet().toArray(
                new String[snpIntervals.size()]);
        Arrays.sort(sortedStrainNames);
        return sortedStrainNames;
    }
    
    /**
     * Convert the intervals to a dataset that JFreeChart can use
     * @param snpIntervals
     *          the intervals
     * @param startInBasePairs
     *          the starting point that we should use
     * @param extentInBasePairs
     *          the extent that we should use
     * @param maximumImageBlockCount
     *          the max number of blocks that we should use for the graph
     * @return
     *          the data set
     */
    private XYDataset snpIntervalsToDataset(
            final Map<String, ? extends List<? extends BasePairInterval>> snpIntervals,
            final long startInBasePairs,
            final long extentInBasePairs,
            final int maximumImageBlockCount)
    {
        DefaultXYZDataset dataset = new DefaultXYZDataset();
        
        final int imageBlockCount;
        final long basePairsPerImageBlock;
        if(maximumImageBlockCount > extentInBasePairs)
        {
            imageBlockCount = (int)extentInBasePairs;
            basePairsPerImageBlock = 1L;
        }
        else
        {
            imageBlockCount = maximumImageBlockCount;
            if(extentInBasePairs % maximumImageBlockCount == 0)
            {
                basePairsPerImageBlock =
                    extentInBasePairs / maximumImageBlockCount;
            }
            else
            {
                basePairsPerImageBlock =
                    (long)Math.ceil((extentInBasePairs / (double)maximumImageBlockCount));
            }
        }
        
        String[] sortedStrainNames = extractSortedStrainNames(snpIntervals);
        
        for(int strainIndex = 0;
            strainIndex < sortedStrainNames.length;
            strainIndex++)
        {
            List<? extends BasePairInterval> snpBlockList =
                snpIntervals.get(sortedStrainNames[strainIndex]);
            int snpBlockIndex = 0;
            BasePairInterval currSnpBlock;
            if(snpBlockList.isEmpty())
            {
                currSnpBlock = null;
            }
            else
            {
                currSnpBlock = snpBlockList.get(snpBlockIndex);
            }
            
            double[][] series = new double[3][imageBlockCount];
            long imageBlockStartInBasePairs = startInBasePairs;
            for(int imageBlockIndex = 0;
                imageBlockIndex < imageBlockCount;
                imageBlockIndex++)
            {
                while(currSnpBlock != null &&
                      currSnpBlock.getEndInBasePairs() < imageBlockStartInBasePairs)
                {
                    snpBlockIndex++;
                    if(snpBlockIndex == snpBlockList.size())
                    {
                        currSnpBlock = null;
                    }
                    else
                    {
                        currSnpBlock = snpBlockList.get(snpBlockIndex);
                    }
                }
                
                final double colorValue;
                if(currSnpBlock == null)
                {
                    // no intersection
                    colorValue = 0.0;
                }
                else
                {
                    if(currSnpBlock.contains(
                            imageBlockStartInBasePairs,
                            basePairsPerImageBlock))
                    {
                        // full containment
                        colorValue = 1.0;
                    }
                    else if(currSnpBlock.intersects(
                            imageBlockStartInBasePairs,
                            basePairsPerImageBlock))
                    {
                        // partial intersection
                        long cumulativeOverlapInBasePairs = 0;
                        for(int innerSnpBlockCursor = snpBlockIndex;
                            innerSnpBlockCursor < snpBlockList.size();
                            innerSnpBlockCursor++)
                        {
                            BasePairInterval innerSnpBlock =
                                snpBlockList.get(innerSnpBlockCursor);
                            long currOverlapInBasePairs =
                                innerSnpBlock.getOverlapInBasePairs(
                                        imageBlockStartInBasePairs,
                                        basePairsPerImageBlock);
                            if(currOverlapInBasePairs == 0)
                            {
                                break;
                            }
                            else
                            {
                                cumulativeOverlapInBasePairs +=
                                    currOverlapInBasePairs;
                            }
                        }
                        
                        double overlapRatio =
                            cumulativeOverlapInBasePairs /
                            (double)basePairsPerImageBlock;
                        colorValue = overlapRatio;
                    }
                    else
                    {
                        // no intersection
                        colorValue = 0.0;
                    }
                }
                
                series[0][imageBlockIndex] = imageBlockStartInBasePairs;
                series[1][imageBlockIndex] = strainIndex;
                series[2][imageBlockIndex] = colorValue;
                
                imageBlockStartInBasePairs += basePairsPerImageBlock;
            }
            dataset.addSeries(sortedStrainNames[strainIndex], series);
        }
        
        return dataset;
    }
}
