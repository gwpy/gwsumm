/*
 * Copyright (C) Duncan Macleod (2013)
 *
 * This file is part of GWSumm
 *
 * GWSumm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GWSumm is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GWSumm.  If not, see <http://www.gnu.org/licenses/>
 */

/* ------------------------------------------------------------------------- */
/* Calendar links                                                            */
var re_dateurl = new RegExp('\/\\d{8}\/')
var re_monthurl = new RegExp('\/\\d{6}\/')

// Move to the date selected
function move_to_date(ev) {
    var newdate = ev.date;
    d = newdate.getDate();
    dd = (d < 10 ? '0' : '') + d;
    m = newdate.getMonth()+1;
    mm = (m < 10 ? '0' : '') + m;
    y = newdate.getFullYear();
    var url = window.location.href;
    if ((ev.viewMode == 'days') && (re_dateurl.test(url))) {
        newurl = url.replace(re_dateurl, '/' + y + mm + dd + '/');
        window.location = newurl;
    }
    else if ((ev.viewMode == 'months') && (re_monthurl.test(url))) {
        newurl = url.replace(re_monthurl, '/' + y + mm + '/');
        window.location = newurl;
    }
    else if ((ev.viewMode == 'years') && (re_yearurl.test(url))) {
        newurl = url.replace(re_yearurl, '/' + y + '/');
        window.location = newurl;
    }
}

// fix width of fixed navbar
function reset_width_on_resize() {
    $('#nav').width($("#nav").width());
}

/* ------------------------------------------------------------------------- */
/* Document ready and loaded                                                  */

// When document is ready, run this stuff:
$(window).load(function() {
    // Set the sticky navigation bar
    $(function() {
        $('#nav-wrapper').height($("#nav").height());
        $('#nav').width($("#nav").width());
        $('#nav').affix({
            offset: { top: $('#nav').offset().top }
        });
    });

    // define the calendar
    $('#calendar').datepicker({weekStart: 1}).on('changeDate', move_to_date);

    // load the fancybox
    $(".fancybox").fancybox({helpers: {title: {type: 'inside'}}});
})

// When document is resized, fix the width of the sticky navbar
var resizeclock;
$(window).resize(function() {
    clearTimeout(resizeclock);
    $('#nav').width($('body').width());
    resizeclock = setTimeout(reset_width_on_resize, 100);
});
