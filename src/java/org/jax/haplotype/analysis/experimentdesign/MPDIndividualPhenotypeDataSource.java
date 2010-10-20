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

package org.jax.haplotype.analysis.experimentdesign;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * A {@link PhenotypeDataSource} that uses a
 * {@link MPDIndividualStrainPhenotypeParser} to get at the data.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class MPDIndividualPhenotypeDataSource implements PhenotypeDataSource
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 1422027862547879921L;

    private static final Logger LOG = Logger.getLogger(
            MPDIndividualPhenotypeDataSource.class.getName());
    
    private static final MPDIndividualStrainPhenotypeParser PHENOTYPE_PARSER =
        new MPDIndividualStrainPhenotypeParser();

    private final File phenotypeFile;
    
    private final String phenotype;
    
    private final Set<String> strainsToAccept;

    private final SexFilter sexToAccept;
    
    private final String name;
    
    /**
     * Constructor
     * @param name
     *          the name of this phenotype data source
     * @param phenotypeFile
     *          the phenotype file (in MPD format)
     * @param phenotype
     *          the phenotype to extract from the file
     * @param strainsToAccept
     *          the strain names to accept
     * @param sexToAccept
     *          the strain sex to accept
     */
    public MPDIndividualPhenotypeDataSource(
            String name,
            File phenotypeFile,
            String phenotype,
            Set<String> strainsToAccept,
            SexFilter sexToAccept)
    {
        this.name = name;
        this.phenotypeFile = phenotypeFile;
        this.phenotype = phenotype;
        this.strainsToAccept = strainsToAccept;
        this.sexToAccept = sexToAccept;
    }
    
    /**
     * {@inheritDoc}
     */
    public String getName()
    {
        return this.name;
    }

    /**
     * {@inheritDoc}
     */
    public Map<String, List<Double>> getPhenotypeData()
    {
        return this.getPhenotypeData(this.strainsToAccept, this.sexToAccept);
    }
    
    /**
     * Get the phenotype using the given constraints
     * @param strainsToAccept
     *          the strains we'll accept
     * @param sexToAccept
     *          the sex we'll accept
     * @return
     *          the strain name to phenotype map
     */
    public Map<String, List<Double>> getPhenotypeData(
            Set<String> strainsToAccept,
            SexFilter sexToAccept)
    {
        try
        {
            return PHENOTYPE_PARSER.parsePhenotypesFromStream(
                    this.phenotype,
                    new FileInputStream(this.phenotypeFile),
                    sexToAccept,
                    strainsToAccept);
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to read phenotype data",
                    ex);
            return null;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return this.getName();
    }
}
