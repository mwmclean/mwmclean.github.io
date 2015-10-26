function(exName, arr){
    $(document).ready(function() {
      setInterval( function() {
	$(arr).map(
	       function( idx , id ){
		   $("#" + id).attr( "src" , exName + "/plots/" + id + ".png" + "?ign=" +
				     (new Date()).valueOf() );
	       } );
	$("#updateMessage").text( "Last refresh was at " + formatted_time() + "." );
      } , 30000 );
    } );
}
