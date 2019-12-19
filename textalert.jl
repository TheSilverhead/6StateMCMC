using PyCall
function TextAlert(number,message)

py"""

def txt(num,mes):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart

    email = "julialangalerts@gmail.com"
    pas = "4"

    sms_gateway = (str(num)+'@tmomail.net')
    smtp = "smtp.gmail.com"
    port = 587
    server = smtplib.SMTP(smtp,port)
    server.starttls()
    server.login(email,pas)

    #MIME module to structure message.
    msg = MIMEMultipart()
    msg['From'] = email
    msg['To'] = sms_gateway
    msg['Subject'] = "Julia Notification\n"
    body = (str(mes)+"\n")
    msg.attach(MIMEText(body, 'plain'))

    sms = msg.as_string()
    server.sendmail(email,sms_gateway,sms)
    server.quit()
"""
py"txt"(number,message)

end
