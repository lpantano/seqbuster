function loadtable(table){

	var srcHolder=document.getElementById('show');
	if ( srcHolder.hasChildNodes() )
	{
		while ( srcHolder.childNodes.length >= 1 )
		{
			srcHolder.removeChild( srcHolder.firstChild );       
		} 
	}
	
	
	var srcTableid=document.getElementById(table);
	
	var srcTable = srcTableid.cloneNode(true);
	srcTable.id = 'visibleTable';
	srcTable.style.visibility='';
	srcHolder.appendChild(srcTable);


}