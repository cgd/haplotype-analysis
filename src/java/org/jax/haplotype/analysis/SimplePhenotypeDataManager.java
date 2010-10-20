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

import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.bind.JAXBContext;

import org.jax.haplotype.analysis.experimentdesign.PhenotypeDataSource;
import org.jax.haplotype.analysis.jaxbfactory.JaxbPhenotypeDataSourceFactory;
import org.jax.haplotype.jaxbgenerated.PhenotypeDataSourceType;
import org.jax.haplotype.jaxbgenerated.PhenotypeDataSources;


/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SimplePhenotypeDataManager
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            SimplePhenotypeDataManager.class.getName());
    
    private static final SimplePhenotypeDataManager instance =
        new SimplePhenotypeDataManager();
    
    private static final String PHENOTYPE_METADATA_RESOURCE_LOCATION =
        "/phenotype-data-sources.xml";
    
    private final Map<String, PhenotypeDataSource> phenotypeDataMap =
        new HashMap<String, PhenotypeDataSource>();
    
    /**
     * 
     */
    private SimplePhenotypeDataManager()
    {
        try
        {
            JAXBContext jaxbContext = JAXBContext.newInstance(
                    PhenotypeDataSources.class.getPackage().getName());
            InputStream phenotypeMetadataStream = SimplePhenotypeDataManager.class.getResourceAsStream(
                    PHENOTYPE_METADATA_RESOURCE_LOCATION);
            PhenotypeDataSources jaxbPhenotypeDataSources =
                (PhenotypeDataSources)jaxbContext.createUnmarshaller().unmarshal(
                        phenotypeMetadataStream);
            
            for(PhenotypeDataSourceType jaxbPhenotypeDataSource:
                jaxbPhenotypeDataSources.getPhenotypeDataSource())
            {
                PhenotypeDataSource phenotypeDataSource =
                    JaxbPhenotypeDataSourceFactory.getPhenotypeDataSource(jaxbPhenotypeDataSource);
                this.phenotypeDataMap.put(
                        phenotypeDataSource.getName(),
                        phenotypeDataSource);
            }
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to read jaxb data",
                    ex);
        }
    }
    
    /**
     * Getter for the instance
     * @return
     *          the instance
     */
    public static SimplePhenotypeDataManager getInstance()
    {
        return SimplePhenotypeDataManager.instance;
    }
    
    /**
     * Getter for the phenotype data map where keys are phenotype names
     * @return
     *          the phenotype data map
     */
    public Map<String, PhenotypeDataSource> getPhenotypeDataMap()
    {
        return this.phenotypeDataMap;
    }
}
