// function showExample(exToShow){
//     window.location.href = "Examples/index.html";
// 	if (document.readyState === "complete") {
// 	    $(exToShow).trigger("click");
// 	}
//     var readyStateCheckInterval = setInterval(function() {
//       if (document.readyState === "complete") {
// 	  clearInterval(readyStateCheckInterval);
// 	  init();
//       }
//     }, 10);
//     	  $(exToShow).trigger("click");
//     // var button = document.getElementById(exToShow);
//     // console.message(button.href);
//     // window.location = button.href;
// }
/* See: http://stackoverflow.com/a/24161846 */
    $(function() {

        // $(document).bind("myCustomEvent", function(e, data) {
        //     alert("Weee, you triggered event id: " + data.EventId + "!");
        // });
	
        // Load the query string specified by the previous page's link
        var exName = getParameterByName('ex');  
        // if(eventId == 'rental') {
        //     // $(document).trigger("myCustomEvent", {EventId: eventId});
	//     $('#rental-btn').trigger('click');
        // } else if (eventId == 'flight'){
	//     $('#'+eventId).trigger('click');
        // } // end if/else
	$('#' + exName + '-btn').trigger('click');
    });

    function getParameterByName(name) {
        name = name.replace(/[\[]/, "\\[").replace(/[\]]/, "\\]");
        var regex = new RegExp("[\\?&]" + name + "=([^&#]*)"),
        results = regex.exec(location.search);
        return results == null ? "" : decodeURIComponent(results[1].replace(/\+/g, " "));
    }
