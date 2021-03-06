/**
 * (C) Copyright IBM Corp. 2010, 2015
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

package com.ibm.bi.dml.hops.rewrite;

import java.util.ArrayList;
import java.util.HashSet;

import com.ibm.bi.dml.hops.DataOp;
import com.ibm.bi.dml.hops.Hop;
import com.ibm.bi.dml.hops.Hop.DataOpTypes;
import com.ibm.bi.dml.hops.HopsException;
import com.ibm.bi.dml.hops.Hop.VisitStatus;

/**
 * This rewrite is a custom rewrite for JMLC in order to replace all persistent reads
 * and writes with transient reads and writes from the symbol table.
 * 
 * 
 */
public class RewriteRemovePersistentReadWrite extends HopRewriteRule
{
	
	private HashSet<String> _inputs = null;
	private HashSet<String> _outputs = null;
	
	public RewriteRemovePersistentReadWrite( String[] in, String[] out )
	{
		_inputs = new HashSet<String>();
		for( String var : in )
			_inputs.add( var );
		_outputs = new HashSet<String>();
		for( String var : out )
			_outputs.add( var );
	}
	
	@Override
	public ArrayList<Hop> rewriteHopDAGs(ArrayList<Hop> roots, ProgramRewriteStatus state)
		throws HopsException
	{
		if( roots == null )
			return null;
		
		for( Hop h : roots ) 
			rule_RemovePersistentDataOp( h );
		
		return roots;
	}

	@Override
	public Hop rewriteHopDAG(Hop root, ProgramRewriteStatus state) 
		throws HopsException
	{
		if( root == null )
			return root;
		
		rule_RemovePersistentDataOp( root );
		
		return root;
	}
	
	/**
	 * 
	 * @param hop
	 * @throws HopsException 
	 */
	private void rule_RemovePersistentDataOp( Hop hop ) 
		throws HopsException
	{
		//check mark processed
		if( hop.getVisited() == VisitStatus.DONE )
			return;
		
		//recursively process childs
		ArrayList<Hop> inputs = hop.getInput();
		for( int i=0; i<inputs.size(); i++ )
			rule_RemovePersistentDataOp( inputs.get(i) );

		//remove cast if unnecessary
		if( hop instanceof DataOp )
		{
			DataOp dop = (DataOp) hop;
			DataOpTypes dotype = dop.getDataOpType();
			
			switch( dotype ) 
			{
				case PERSISTENTREAD:
					if( _inputs.contains(dop.getName()) )
						dop.setDataOpType(DataOpTypes.TRANSIENTREAD);
					break;
				case PERSISTENTWRITE:
					if( _outputs.contains(dop.getName()) )
						dop.setDataOpType(DataOpTypes.TRANSIENTWRITE);
					break;
				default:
					//do nothing
			}
		}
		
		//mark processed
		hop.setVisited( VisitStatus.DONE );
	}
}
