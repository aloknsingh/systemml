package com.ibm.bi.dml.utils;

@SuppressWarnings("serial")
public class LopsException extends DMLException 
{
	
	public LopsException(String message)
	{
		super(message);
	}

	public LopsException(Exception e) {
		super(e);
	}
	
	public LopsException(String message, Throwable cause) {
	    super(message, cause);
	}

}
