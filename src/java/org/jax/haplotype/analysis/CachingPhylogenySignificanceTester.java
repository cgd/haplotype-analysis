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

package org.jax.haplotype.analysis;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.inference.CachingGenomeDataManager;
import org.jax.haplotype.inference.CachingPhylogenyDataManager;
import org.jax.haplotype.phylogeny.data.NoValidPhylogenyException;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTestResult;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdgeWithRealValue;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CachingPhylogenySignificanceTester
{
    private static final Logger LOG = Logger.getLogger(
            CachingPhylogenySignificanceTester.class.getName());
    
    private final PhylogenySignificanceTester phylogenySignificanceTester =
        new PhylogenySignificanceTester();
    
    private final SimplePhenotypeDataManager phenotypeDataManager;
    
    private final CachingPhylogenyDataManager phylogenyDataManager;
    
    private final Map<String, File> resultsCacheFileMap =
        new HashMap<String, File>();

    private static final CachingPhylogenySignificanceTester instance =
        new CachingPhylogenySignificanceTester(
                CachingPhylogenyDataManager.getInstance(),
                SimplePhenotypeDataManager.getInstance());
    
    protected static final String CONCATINATION_STRING = "&";
    
    /**
     * Get the singleton instance
     * @return
     *          the instance
     */
    public static CachingPhylogenySignificanceTester getInstance()
    {
        return CachingPhylogenySignificanceTester.instance;
    }
    
    /**
     * @param phenotypeDataManager
     * @param phylogenyDataManager
     */
    private CachingPhylogenySignificanceTester(
            CachingPhylogenyDataManager phylogenyDataManager,
            SimplePhenotypeDataManager phenotypeDataManager)
    {
        this.phenotypeDataManager = phenotypeDataManager;
        this.phylogenyDataManager = phylogenyDataManager;
    }
    
    /**
     * Get the phylogeny test results
     * @param phenotypeName
     *          the phenotype name
     * @param genomeName
     *          the genotype name
     * @param strainNames
     *          the strain names
     * @param chromosomeNumber
     *          the chromosome number
     * @return
     *          the test results
     * @throws NoValidPhylogenyException 
     * @throws IOException 
     */
    @SuppressWarnings("unchecked")
    public synchronized List<PhylogenyTestResult> getPhylogenyTestResults(
            String phenotypeName,
            String genomeName,
            String[] strainNames,
            int chromosomeNumber) throws IOException, NoValidPhylogenyException
    {
        File cacheFile = this.getCacheFile(
                phenotypeName,
                genomeName,
                strainNames,
                chromosomeNumber);
        
        if(cacheFile.createNewFile())
        {
            cacheFile.deleteOnExit();
            
            List<PhylogenyInterval> phylogenyIntervals = this.phylogenyDataManager.getPhylogeneticIntervals(
                    genomeName,
                    strainNames,
                    chromosomeNumber);
            PhenotypeDataSource phenotypeDataSource = this.phenotypeDataManager.getPhenotypeDataMap().get(
                    phenotypeName);
            
            Map<String, List<Double>> phenotypeData =
                phenotypeDataSource.getPhenotypeData();
            phenotypeData.keySet().retainAll(Arrays.asList(strainNames));
            CachingGenomeDataManager genomeDataManager =
                this.phylogenyDataManager.getGenomeDataManager();
            GenomeDataSource selectedGenome =
                genomeDataManager.getGenomeDataMap().get(genomeName);
            phenotypeData.keySet().retainAll(selectedGenome.getAvailableStrains());
            
            List<PhylogenyTestResult> testResults = new ArrayList<PhylogenyTestResult>(
                    phylogenyIntervals.size());
            
            for(PhylogenyInterval phylogenyInterval: phylogenyIntervals)
            {
                PhylogenyTreeNode prunedPhylogeny = phylogenyInterval.getPhylogeny().createStrainPrunedTree(
                        phenotypeData.keySet());
                PhylogenyTreeNode treeWithSignificance =
                    this.phylogenySignificanceTester.testMultipleResponseSignificance(
                            prunedPhylogeny,
                            phenotypeData);
                PhylogenyTreeEdgeWithRealValue minPValueEdge =
                    PhylogenyTreeEdgeWithRealValue.getEdgeWithMininumValue(
                            treeWithSignificance);
                PhylogenyTestResult phylogenyTestResult = new PhylogenyTestResult(
                        new PhylogenyInterval(prunedPhylogeny, phylogenyInterval.getInterval()),
                        minPValueEdge == null ? 1.0 : minPValueEdge.getRealValue());
                testResults.add(phylogenyTestResult);
            }
            
            ObjectOutputStream objectOut = new ObjectOutputStream(
                    new BufferedOutputStream(new FileOutputStream(cacheFile)));
            objectOut.writeObject(testResults);
            objectOut.flush();
            objectOut.close();
            
            return testResults;
        }
        else
        {
            // the file is already cached. load it and return it
            ObjectInputStream objectIn = new ObjectInputStream(
                    new BufferedInputStream(new FileInputStream(cacheFile)));
            
            try
            {
                List<PhylogenyTestResult> testResults =
                    (List<PhylogenyTestResult>)objectIn.readObject();
                return testResults;
            }
            catch(ClassNotFoundException ex)
            {
                LOG.log(Level.SEVERE,
                        "failed to read in cached phylogeny results",
                        ex);
                return null;
            }
        }
    }

    /**
     * Get a cache file to be used for the given parameters. This function
     * will not create the file on disk
     */
    private File getCacheFile(
            String phenotypeName,
            String genomeName,
            String[] strainNames,
            int chromosomeNumber)
    {
        strainNames = strainNames.clone();
        Arrays.sort(strainNames);
        
        String cacheKeyString =
            phenotypeName + CONCATINATION_STRING +
            genomeName + CONCATINATION_STRING +
            Arrays.toString(strainNames) + CONCATINATION_STRING +
            chromosomeNumber;
        
        File cacheFile = this.resultsCacheFileMap.get(cacheKeyString);
        if(cacheFile == null)
        {
            File directory = new File(System.getProperty("java.io.tmpdir"));
            cacheFile = new File(
                    directory,
                    "bham-cache-" + this.resultsCacheFileMap.size() + ".bahm");
            this.resultsCacheFileMap.put(
                    cacheKeyString,
                    cacheFile);
        }
        return cacheFile;
    }
}
