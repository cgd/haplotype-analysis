/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.analysis.experimentdesign;

import java.io.Serializable;
import java.util.List;
import java.util.Map;

/**
 * An interface definition for a phenotype data source
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface PhenotypeDataSource extends Serializable
{
    /**
     * Get the phenotyp data
     * @return
     *          the strain name to phenotype map
     */
    public Map<String, List<Double>> getPhenotypeData();

    /**
     * Get the name of this phenotype data source
     * @return
     *          the name
     */
    public String getName();
}
