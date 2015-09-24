// Miscellaneous javascript functions.

// Shows the hours:mins:secs with AM or PM
  function formatted_time() {
    var currentTime = new Date()
    var hours = currentTime.getHours()
    var minutes = currentTime.getMinutes()
    var seconds = currentTime.getSeconds()

    var suffix = "AM";
    if (hours >= 12) {
        suffix = "PM";
        hours = hours - 12;
    }

    if (hours == 0) {
        hours = 12;
    }

    if (minutes < 10)
        minutes = "0" + minutes

    if (seconds < 10)
        seconds = "0" + seconds 

    return hours + ":" + minutes + ":" + seconds + " " + suffix
  }


