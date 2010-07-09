/**
 * General class to encapsulate exceptions having to do with MaeNodes, such as
 * some method being invoked on the wrong part of the .mae node tree.
 * 
 * @author Tom Joseph <thomas.joseph@mssm.edu>
 */

public class MaeNodeException extends Exception
{
	public MaeNodeException(String string)
	{
		super(string);
	}

}
